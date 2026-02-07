/* BEGIN Header */
/**
 ******************************************************************************
 * \file            AHRS_VQF.c
 * \author          Andrea Vivani
 * \brief           Implementation of VQF for attitude and heading estimation
 ******************************************************************************
 * \copyright
 *
 * Copyright 2026 Andrea Vivani
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the “Software”), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 *
 ******************************************************************************
 */
/* END Header */

/* Configuration check -------------------------------------------------------*/
#if !defined(ADVUTILS_USE_DYNAMIC_ALLOCATION)
#error ADVUTILS_USE_DYNAMIC_ALLOCATION must be set to use AHRS_VQF (dynamic matrix allocation)
#endif

/* Includes ------------------------------------------------------------------*/

#include "AHRS_VQF.h"

#include <float.h>
#include <math.h>
#include <string.h>

#include "basicMath.h"

/* Macros --------------------------------------------------------------------*/

#ifndef VQF_EPS
#define VQF_EPS (FLT_EPSILON)
#endif

/* Typedefs ------------------------------------------------------------------*/

typedef struct {
    /* Core filter time constants [s] */
    float tauAcc;
    float tauMag;

    /* Feature toggles */
    uint8_t motionBiasEstEnabled;
    uint8_t restBiasEstEnabled;
    uint8_t magDistRejectionEnabled;

    /* Gyro bias estimation tuning */
    float biasSigmaInit;
    float biasForgettingTime;
    float biasClip;
    float biasSigmaMotion;
    float biasVerticalForgettingFactor;
    float biasSigmaRest;

    /* Rest detection */
    float restMinT;
    float restFilterTau;
    float restThGyr;
    float restThAcc;

    /* Magnetic disturbance rejection */
    float magCurrentTau;
    float magRefTau;
    float magNormTh;
    float magDipTh;
    float magNewTime;
    float magNewFirstTime;
    float magNewMinGyr;
    float magMinUndisturbedTime;
    float magMaxRejectionTime;
    float magRejectionFactor;
} AHRS_VQF_Params_t;

typedef struct {
    float gyrTs;
    float accTs;
    float magTs;

    matrix_t accLpB;
    matrix_t accLpA;

    float kMag;

    float biasP0;
    float biasV;
    float biasMotionW;
    float biasVerticalW;
    float biasRestW;

    matrix_t restGyrLpB;
    matrix_t restGyrLpA;
    matrix_t restAccLpB;
    matrix_t restAccLpA;

    float kMagRef;
    matrix_t magNormDipLpB;
    matrix_t magNormDipLpA;
} AHRS_VQF_Coeffs_t;

typedef struct {
    quaternion_t gyrQuat;
    quaternion_t accQuat;
    float delta;

    uint8_t restDetected;
    uint8_t magDistDetected;

    axis3f_t lastAccLp;
    matrix_t accLpState;
    float lastAccCorrAngularRate;

    float kMagInit;
    float lastMagDisAngle;
    float lastMagCorrAngularRate;

    axis3f_t bias;
    matrix_t biasP;

    matrix_t motionBiasEstRLpState;
    matrix_t motionBiasEstBiasLpState;

    matrix_t restLastSquaredDeviations;
    float restT;
    axis3f_t restLastGyrLp;
    matrix_t restGyrLpState;
    axis3f_t restLastAccLp;
    matrix_t restAccLpState;

    float magRefNorm;
    float magRefDip;
    float magUndisturbedT;
    float magRejectT;

    float magCandidateNorm;
    float magCandidateDip;
    float magCandidateT;

    matrix_t magNormDip;
    matrix_t magNormDipLpState;
} AHRS_VQF_State_t;

/* Private variables ---------------------------------------------------------*/

AHRS_VQF_Params_t params;
AHRS_VQF_Coeffs_t coeffs;
AHRS_VQF_State_t state;

/* Scratch matrices (allocated once in init) */
matrix_t _R;
matrix_t _TMP33a;
matrix_t _TMP33b;
matrix_t _TMP33c;
matrix_t _TMP33d;
matrix_t _TMP31a;
matrix_t _TMP31b;
matrix_t _TMP21a;
matrix_t _e;

/* Functions -----------------------------------------------------------------*/

static inline float vqf_square(float x) { return x * x; }

static inline float vqf_wrapPi(float angle) {
    if (angle > constPI) {
        angle -= 2.0f * constPI;
    } else if (angle < -constPI) {
        angle += 2.0f * constPI;
    }
    return angle;
}

static inline axis3f_t vqf_bodyNedToEnu(axis3f_t v_ned) {
    axis3f_t v_enu;
    v_enu.x = v_ned.y;
    v_enu.y = v_ned.x;
    v_enu.z = -v_ned.z;
    return v_enu;
}

static inline axis3f_t vqf_bodyEnuToNed(axis3f_t v_enu) {
    axis3f_t v_ned;
    v_ned.x = v_enu.y;
    v_ned.y = v_enu.x;
    v_ned.z = -v_enu.z;
    return v_ned;
}

static void vqf_matrixFill(matrix_t* m, float value) {
    if (!m || !m->data) {
        return;
    }
    for (uint16_t i = 0; i < (uint16_t)m->rows * (uint16_t)m->cols; i++) {
        m->data[i] = value;
    }
}

static float vqf_gainFromTau(float tau, float Ts) {
    if (tau < 0.0f) {
        return 0.0f;
    }
    if (tau == 0.0f) {
        return 1.0f;
    }
    return 1.0f - expf(-Ts / tau);
}

static void vqf_filterCoeffs(float tau, float Ts, matrix_t* outB_3x1, matrix_t* outA_2x1) {
    /* tau <= 0 means "no filtering" (pass-through). This is consistent with the
     * behaviour of gainFromTau and avoids numerical issues for very small tau.
     */
    if (tau <= 0.0f) {
        outB_3x1->data[0] = 1.0f;
        outB_3x1->data[1] = 0.0f;
        outB_3x1->data[2] = 0.0f;
        outA_2x1->data[0] = 0.0f;
        outA_2x1->data[1] = 0.0f;
        return;
    }

    /* 2nd order Butterworth LPF based on https://stackoverflow.com/a/52764064 */
    const float sqrt2 = SQRT(2.0f);
    const float fc = (sqrt2 / (2.0f * constPI)) / tau;
    const float C = TAN(constPI * fc * Ts);
    const float D = C * C + sqrt2 * C + 1.0f;

    const float b0 = (C * C) / D;
    outB_3x1->data[0] = b0;
    outB_3x1->data[1] = 2.0f * b0;
    outB_3x1->data[2] = b0;

    outA_2x1->data[0] = 2.0f * (C * C - 1.0f) / D;      /* a1 (a0 = 1) */
    outA_2x1->data[1] = (1.0f - sqrt2 * C + C * C) / D; /* a2 */
}

static void vqf_filterInitialState(float x0, const matrix_t* b_3x1, const matrix_t* a_2x1, float* out2) {
    out2[0] = x0 * (1.0f - b_3x1->data[0]);
    out2[1] = x0 * (b_3x1->data[2] - a_2x1->data[1]);
}

static float vqf_filterStep(float x, const matrix_t* b_3x1, const matrix_t* a_2x1, float* state2) {
    /* Difference equations based on scipy.signal.lfilter documentation (a0 == 1) */
    const float y = b_3x1->data[0] * x + state2[0];
    state2[0] = b_3x1->data[1] * x - a_2x1->data[0] * y + state2[1];
    state2[1] = b_3x1->data[2] * x - a_2x1->data[1] * y;
    return y;
}

static void vqf_filterVec(const matrix_t* x_Nx1, float tau, float Ts, const matrix_t* b_3x1, const matrix_t* a_2x1, matrix_t* state_Nx2, matrix_t* out_Nx1) {
    /* state_Nx2 is also used as initialization scratch, following the original VQF implementation */
    const uint16_t N = (uint16_t)x_Nx1->rows * (uint16_t)x_Nx1->cols;
    float* st = state_Nx2->data;
    float* out = out_Nx1->data;
    const float* x = x_Nx1->data;

    if (isnan(st[0])) {
        /* Initialization phase: average the first samples for duration tau */
        if (isnan(st[1])) {
            st[1] = 0.0f; /* sample count */
            for (uint16_t i = 0; i < N; i++) {
                st[2u + i] = 0.0f; /* sums */
            }
        }

        st[1] += 1.0f;
        for (uint16_t i = 0; i < N; i++) {
            st[2u + i] += x[i];
            out[i] = st[2u + i] / st[1];
        }

        if (st[1] * Ts >= tau) {
            for (uint16_t i = 0; i < N; i++) {
                vqf_filterInitialState(out[i], b_3x1, a_2x1, &st[2u * i]);
            }
        }
        return;
    }

    for (uint16_t i = 0; i < N; i++) {
        out[i] = vqf_filterStep(x[i], b_3x1, a_2x1, &st[2u * i]);
    }
}

static axis3f_t vqf_vecNormalize(axis3f_t v) {
    const float n2 = v.x * v.x + v.y * v.y + v.z * v.z;
    if (n2 <= VQF_EPS) {
        return v;
    }
    const float inv_n = INVSQRT(n2);
    v.x *= inv_n;
    v.y *= inv_n;
    v.z *= inv_n;
    return v;
}

static float vqf_vecNorm(axis3f_t v) { return SQRT(v.x * v.x + v.y * v.y + v.z * v.z); }

static axis3f_t vqf_quatRotateForward(const quaternion_t* q, axis3f_t v) {
    /* quaternionRotation() implements q_conj * v * q. To get q * v * q_conj, rotate with q_conj. */
    quaternion_t qc;
    quaternionConj((quaternion_t*)q, &qc);

    quaternion_t qv;
    quaternion_t qo;
    qv.q0 = 0.0f;
    qv.q1 = v.x;
    qv.q2 = v.y;
    qv.q3 = v.z;
    quaternionRotation(&qc, &qv, &qo);

    axis3f_t out;
    out.x = qo.q1;
    out.y = qo.q2;
    out.z = qo.q3;
    return out;
}

static void vqf_applyDelta(quaternion_t* q, float delta_rad) {
    /* q <- q_delta ⊗ q, q_delta = [cos(d/2), 0, 0, sin(d/2)] */
    quaternion_t qd;
    const float half = 0.5f * delta_rad;
    qd.q0 = COS(half);
    qd.q1 = 0.0f;
    qd.q2 = 0.0f;
    qd.q3 = SIN(half);

    quaternion_t tmp;
    quaternionMult(&qd, q, &tmp);
    *q = tmp;
}

static void vqf_getQuat3D_enu(quaternion_t* q_out) { *q_out = state.gyrQuat; }

static void vqf_getQuat6D_enu(quaternion_t* q_out) { quaternionMult((quaternion_t*)&state.accQuat, (quaternion_t*)&state.gyrQuat, q_out); }

static void vqf_getQuat9D_enu(quaternion_t* q_out) {
    vqf_getQuat6D_enu(q_out);
    vqf_applyDelta(q_out, state.delta);
    quaternionNorm(q_out);
}

static void vqf_quatEnuToNed(const quaternion_t* q_enu, quaternion_t* q_ned) {
    /* qP is 180° rotation about axis (1,1,0)/sqrt(2). Its own inverse. */
    const float inv_sqrt2 = INVSQRT(2.0f);
    quaternion_t qP;
    qP.q0 = 0.0f;
    qP.q1 = inv_sqrt2;
    qP.q2 = inv_sqrt2;
    qP.q3 = 0.0f;

    quaternion_t tmp;
    quaternionMult((quaternion_t*)&qP, (quaternion_t*)q_enu, &tmp);
    quaternionMult(&tmp, (quaternion_t*)&qP, q_ned);
    quaternionNorm(q_ned);
}

static void vqf_setMatIdentityScaled(matrix_t* m3x3, float scale) {
    matrixZeros(m3x3);
    ELEMP(m3x3, 0, 0) = scale;
    ELEMP(m3x3, 1, 1) = scale;
    ELEMP(m3x3, 2, 2) = scale;
}

static void vqf_initDefaultParams(AHRS_VQF_Params_t* p) {
    p->tauAcc = 3.0f;
    p->tauMag = 9.0f;

    p->motionBiasEstEnabled = 1u;
    p->restBiasEstEnabled = 1u;
    p->magDistRejectionEnabled = 1u;

    p->biasSigmaInit = 0.5f;
    p->biasForgettingTime = 100.0f;
    p->biasClip = 2.0f;
    p->biasSigmaMotion = 0.1f;
    p->biasVerticalForgettingFactor = 0.0001f;
    p->biasSigmaRest = 0.03f;

    p->restMinT = 1.5f;
    p->restFilterTau = 0.5f;
    p->restThGyr = 2.0f;
    p->restThAcc = 0.5f;

    p->magCurrentTau = 0.05f;
    p->magRefTau = 20.0f;
    p->magNormTh = 0.1f;
    p->magDipTh = 10.0f;
    p->magNewTime = 20.0f;
    p->magNewFirstTime = 5.0f;
    p->magNewMinGyr = 20.0f;
    p->magMinUndisturbedTime = 0.5f;
    p->magMaxRejectionTime = 60.0f;
    p->magRejectionFactor = 2.0f;
}

static void vqf_updateCoeffs() {
    AHRS_VQF_Params_t* p = &params;
    AHRS_VQF_Coeffs_t* c = &coeffs;

    vqf_filterCoeffs(p->tauAcc, c->accTs, &c->accLpB, &c->accLpA);
    c->kMag = vqf_gainFromTau(p->tauMag, c->magTs);

    /* Bias covariance is kept in units of (0.01 deg/s)^2 (same as original VQF) */
    c->biasP0 = vqf_square(p->biasSigmaInit * 100.0f);
    c->biasV = vqf_square(0.1f * 100.0f) * c->accTs / p->biasForgettingTime;

    const float pMotion = vqf_square(p->biasSigmaMotion * 100.0f);
    c->biasMotionW = vqf_square(pMotion) / c->biasV + pMotion;
    c->biasVerticalW = c->biasMotionW / fmaxf(p->biasVerticalForgettingFactor, 1e-10f);

    const float pRest = vqf_square(p->biasSigmaRest * 100.0f);
    c->biasRestW = vqf_square(pRest) / c->biasV + pRest;

    vqf_filterCoeffs(p->restFilterTau, c->gyrTs, &c->restGyrLpB, &c->restGyrLpA);
    vqf_filterCoeffs(p->restFilterTau, c->accTs, &c->restAccLpB, &c->restAccLpA);

    c->kMagRef = vqf_gainFromTau(p->magRefTau, c->magTs);
    if (p->magCurrentTau > 0.0f) {
        vqf_filterCoeffs(p->magCurrentTau, c->magTs, &c->magNormDipLpB, &c->magNormDipLpA);
    } else {
        vqf_matrixFill(&c->magNormDipLpB, NAN);
        vqf_matrixFill(&c->magNormDipLpA, NAN);
    }
}

static void vqf_resetState() {
    /* Quaternions */
    state.gyrQuat.q0 = 1.0f;
    state.gyrQuat.q1 = 0.0f;
    state.gyrQuat.q2 = 0.0f;
    state.gyrQuat.q3 = 0.0f;

    state.accQuat.q0 = 1.0f;
    state.accQuat.q1 = 0.0f;
    state.accQuat.q2 = 0.0f;
    state.accQuat.q3 = 0.0f;

    state.delta = 0.0f;

    state.restDetected = 0u;
    state.magDistDetected = 1u;

    state.lastAccLp.x = 0.0f;
    state.lastAccLp.y = 0.0f;
    state.lastAccLp.z = 0.0f;
    vqf_matrixFill(&state.accLpState, NAN);
    state.lastAccCorrAngularRate = 0.0f;

    state.kMagInit = 1.0f;
    state.lastMagDisAngle = 0.0f;
    state.lastMagCorrAngularRate = 0.0f;

    state.bias.x = 0.0f;
    state.bias.y = 0.0f;
    state.bias.z = 0.0f;
    vqf_setMatIdentityScaled(&state.biasP, coeffs.biasP0);

    /* Motion bias estimator LPFs */
    vqf_matrixFill(&state.motionBiasEstRLpState, NAN);
    vqf_matrixFill(&state.motionBiasEstBiasLpState, NAN);

    /* Rest detection */
    vqf_matrixFill(&state.restLastSquaredDeviations, 0.0f);
    state.restT = 0.0f;
    state.restLastGyrLp.x = 0.0f;
    state.restLastGyrLp.y = 0.0f;
    state.restLastGyrLp.z = 0.0f;
    vqf_matrixFill(&state.restGyrLpState, NAN);
    state.restLastAccLp.x = 0.0f;
    state.restLastAccLp.y = 0.0f;
    state.restLastAccLp.z = 0.0f;
    vqf_matrixFill(&state.restAccLpState, NAN);

    /* Magnetic reference and disturbance rejection */
    state.magRefNorm = 0.0f;
    state.magRefDip = 0.0f;
    state.magUndisturbedT = 0.0f;
    state.magRejectT = params.magMaxRejectionTime;

    state.magCandidateNorm = -1.0f;
    state.magCandidateDip = 0.0f;
    state.magCandidateT = 0.0f;

    vqf_matrixFill(&state.magNormDip, 0.0f);
    vqf_matrixFill(&state.magNormDipLpState, NAN);
}

/*------------------------------------Public functions--------------------------------------------*/

/*------------------------------------Initialization--------------------------------------------*/

void AHRS_VQF_Init(float gyrTs_s, float accTs_s, float magTs_s) {

    vqf_initDefaultParams(&params);

    /* Sampling times */
    coeffs.gyrTs = gyrTs_s;
    coeffs.accTs = (accTs_s > 0.0f) ? accTs_s : gyrTs_s;
    coeffs.magTs = (magTs_s > 0.0f) ? magTs_s : gyrTs_s;

    /* Allocate coefficient matrices */
    matrixInit(&coeffs.accLpB, 3, 1);
    matrixInit(&coeffs.accLpA, 2, 1);
    matrixInit(&coeffs.restGyrLpB, 3, 1);
    matrixInit(&coeffs.restGyrLpA, 2, 1);
    matrixInit(&coeffs.restAccLpB, 3, 1);
    matrixInit(&coeffs.restAccLpA, 2, 1);
    matrixInit(&coeffs.magNormDipLpB, 3, 1);
    matrixInit(&coeffs.magNormDipLpA, 2, 1);

    /* Allocate state matrices */
    matrixInit(&state.accLpState, 3, 2);
    matrixInit(&state.biasP, 3, 3);
    matrixInit(&state.motionBiasEstRLpState, 9, 2);
    matrixInit(&state.motionBiasEstBiasLpState, 2, 2);
    matrixInit(&state.restLastSquaredDeviations, 2, 1);
    matrixInit(&state.restGyrLpState, 3, 2);
    matrixInit(&state.restAccLpState, 3, 2);
    matrixInit(&state.magNormDip, 2, 1);
    matrixInit(&state.magNormDipLpState, 2, 2);

    /* Allocate scratch matrices */
    matrixInit(&_R, 3, 3);
    matrixInit(&_TMP33a, 3, 3);
    matrixInit(&_TMP33b, 3, 3);
    matrixInit(&_TMP33c, 3, 3);
    matrixInit(&_TMP33d, 3, 3);
    matrixInit(&_TMP31a, 3, 1);
    matrixInit(&_TMP31b, 3, 1);
    matrixInit(&_TMP21a, 2, 1);
    matrixInit(&_e, 3, 1);

    vqf_updateCoeffs();
    vqf_resetState();
}

void AHRS_VQF_Deinit() {

    /* Coeff matrices */
    matrixDelete(&coeffs.accLpB);
    matrixDelete(&coeffs.accLpA);
    matrixDelete(&coeffs.restGyrLpB);
    matrixDelete(&coeffs.restGyrLpA);
    matrixDelete(&coeffs.restAccLpB);
    matrixDelete(&coeffs.restAccLpA);
    matrixDelete(&coeffs.magNormDipLpB);
    matrixDelete(&coeffs.magNormDipLpA);

    /* State matrices */
    matrixDelete(&state.accLpState);
    matrixDelete(&state.biasP);
    matrixDelete(&state.motionBiasEstRLpState);
    matrixDelete(&state.motionBiasEstBiasLpState);
    matrixDelete(&state.restLastSquaredDeviations);
    matrixDelete(&state.restGyrLpState);
    matrixDelete(&state.restAccLpState);
    matrixDelete(&state.magNormDip);
    matrixDelete(&state.magNormDipLpState);

    /* Scratch matrices */
    matrixDelete(&_R);
    matrixDelete(&_TMP33a);
    matrixDelete(&_TMP33b);
    matrixDelete(&_TMP33c);
    matrixDelete(&_TMP33d);
    matrixDelete(&_TMP31a);
    matrixDelete(&_TMP31b);
    matrixDelete(&_TMP21a);
    matrixDelete(&_e);
}

void AHRS_VQF_Reset() {
    vqf_updateCoeffs();
    vqf_resetState();
}

/*------------------------------------Update functions--------------------------------------------*/

void AHRS_VQF_UpdateGyro(axis3f_t gyr_rad_s) {

    /* Convert from body-NED to internal body-ENU */
    const axis3f_t gyr = vqf_bodyNedToEnu(gyr_rad_s);

    /* Rest detection based on gyro */
    if (params.restBiasEstEnabled || params.magDistRejectionEnabled) {
        matrix_t x = {(float*)&gyr, 3, 1};
        matrix_t out = {(float*)&state.restLastGyrLp, 3, 1};
        vqf_filterVec(&x, params.restFilterTau, coeffs.gyrTs, &coeffs.restGyrLpB, &coeffs.restGyrLpA, &state.restGyrLpState, &out);

        const float dx = gyr.x - state.restLastGyrLp.x;
        const float dy = gyr.y - state.restLastGyrLp.y;
        const float dz = gyr.z - state.restLastGyrLp.z;
        state.restLastSquaredDeviations.data[0] = dx * dx + dy * dy + dz * dz;

        const float biasClipRad = params.biasClip * (constPI / 180.0f);
        const float thGyrRad = params.restThGyr * (constPI / 180.0f);
        if (state.restLastSquaredDeviations.data[0] >= thGyrRad * thGyrRad || fabsf(state.restLastGyrLp.x) > biasClipRad
            || fabsf(state.restLastGyrLp.y) > biasClipRad || fabsf(state.restLastGyrLp.z) > biasClipRad) {
            state.restT = 0.0f;
            state.restDetected = 0u;
        }
    }

    /* Remove estimated bias (internal ENU) */
    axis3f_t gyrNoBias;
    gyrNoBias.x = gyr.x - state.bias.x;
    gyrNoBias.y = gyr.y - state.bias.y;
    gyrNoBias.z = gyr.z - state.bias.z;

    /* Gyro integration */
    const float gyrNorm = vqf_vecNorm(gyrNoBias);
    const float angle = gyrNorm * coeffs.gyrTs;
    if (gyrNorm > VQF_EPS) {
        quaternion_t dq;
        const float half = 0.5f * angle;
        const float s = SIN(half) / gyrNorm;
        dq.q0 = COS(half);
        dq.q1 = s * gyrNoBias.x;
        dq.q2 = s * gyrNoBias.y;
        dq.q3 = s * gyrNoBias.z;
        quaternionMult(&state.gyrQuat, &dq, &state.gyrQuat);
        quaternionNorm(&state.gyrQuat);
    }
}

void AHRS_VQF_UpdateAcc(axis3f_t acc) {

    /* Convert from body-NED to internal body-ENU */
    const axis3f_t accBody = vqf_bodyNedToEnu(acc);

    /* Ignore [0 0 0] */
    if (accBody.x == 0.0f && accBody.y == 0.0f && accBody.z == 0.0f) {
        return;
    }

    /* Rest detection based on accelerometer */
    if (params.restBiasEstEnabled) {
        matrix_t x = {(float*)&accBody, 3, 1};
        matrix_t out = {(float*)&state.restLastAccLp, 3, 1};
        vqf_filterVec(&x, params.restFilterTau, coeffs.accTs, &coeffs.restAccLpB, &coeffs.restAccLpA, &state.restAccLpState, &out);

        const float dx = accBody.x - state.restLastAccLp.x;
        const float dy = accBody.y - state.restLastAccLp.y;
        const float dz = accBody.z - state.restLastAccLp.z;
        state.restLastSquaredDeviations.data[1] = dx * dx + dy * dy + dz * dz;

        if (state.restLastSquaredDeviations.data[1] >= params.restThAcc * params.restThAcc) {
            state.restT = 0.0f;
            state.restDetected = 0u;
        } else {
            state.restT += coeffs.accTs;
            if (state.restT >= params.restMinT) {
                state.restDetected = 1u;
            }
        }
    }

    /* Filter acc in inertial frame (internal ENU) */
    axis3f_t accEarth = vqf_quatRotateForward(&state.gyrQuat, accBody);

    matrix_t accEarthMat = {(float*)&accEarth, 3, 1};
    matrix_t lastAccLpMat = {(float*)&state.lastAccLp, 3, 1};
    vqf_filterVec(&accEarthMat, params.tauAcc, coeffs.accTs, &coeffs.accLpB, &coeffs.accLpA, &state.accLpState, &lastAccLpMat);

    /* Transform to 6D earth frame and normalize */
    accEarth = vqf_quatRotateForward(&state.accQuat, state.lastAccLp);
    accEarth = vqf_vecNormalize(accEarth);

    /* Inclination correction */
    quaternion_t accCorr;
    const float qw = SQRT((accEarth.z + 1.0f) * 0.5f);
    if (qw > 1e-6f) {
        accCorr.q0 = qw;
        accCorr.q1 = 0.5f * accEarth.y / qw;
        accCorr.q2 = -0.5f * accEarth.x / qw;
        accCorr.q3 = 0.0f;
    } else {
        /* Avoid numeric issues when acc is close to [0 0 -1] */
        accCorr.q0 = 0.0f;
        accCorr.q1 = 1.0f;
        accCorr.q2 = 0.0f;
        accCorr.q3 = 0.0f;
    }

    quaternionMult(&accCorr, &state.accQuat, &state.accQuat);
    quaternionNorm(&state.accQuat);

    state.lastAccCorrAngularRate = acosf(CONSTRAIN(accEarth.z, -1.0f, 1.0f)) / coeffs.accTs;

    /* Gyro bias estimation */
    if (params.motionBiasEstEnabled || params.restBiasEstEnabled) {
        const float biasClipRad = params.biasClip * (constPI / 180.0f);

        quaternion_t accGyrQuat;
        vqf_getQuat6D_enu(&accGyrQuat);

        /* Build rotation matrix R(accGyrQuat) in internal ENU */
        const float q0 = accGyrQuat.q0;
        const float q1 = accGyrQuat.q1;
        const float q2 = accGyrQuat.q2;
        const float q3 = accGyrQuat.q3;

        ELEM(_R, 0, 0) = 1.0f - 2.0f * (q2 * q2 + q3 * q3);
        ELEM(_R, 0, 1) = 2.0f * (q2 * q1 - q0 * q3);
        ELEM(_R, 0, 2) = 2.0f * (q0 * q2 + q3 * q1);
        ELEM(_R, 1, 0) = 2.0f * (q0 * q3 + q2 * q1);
        ELEM(_R, 1, 1) = 1.0f - 2.0f * (q1 * q1 + q3 * q3);
        ELEM(_R, 1, 2) = 2.0f * (q2 * q3 - q1 * q0);
        ELEM(_R, 2, 0) = 2.0f * (q3 * q1 - q0 * q2);
        ELEM(_R, 2, 1) = 2.0f * (q0 * q1 + q3 * q2);
        ELEM(_R, 2, 2) = 1.0f - 2.0f * (q1 * q1 + q2 * q2);

        /* biasLp = (R * bias)_{x,y} */
        matrix_t biasVec = {(float*)&state.bias, 3, 1};
        matrixMult(&_R, &biasVec, &_TMP31a);
        _TMP21a.data[0] = _TMP31a.data[0];
        _TMP21a.data[1] = _TMP31a.data[1];

        /* Low-pass filter R and biasLp */
        matrix_t Rvec = {_R.data, 9, 1};
        vqf_filterVec(&Rvec, params.tauAcc, coeffs.accTs, &coeffs.accLpB, &coeffs.accLpA, &state.motionBiasEstRLpState, &Rvec);
        vqf_filterVec(&_TMP21a, params.tauAcc, coeffs.accTs, &coeffs.accLpB, &coeffs.accLpA, &state.motionBiasEstBiasLpState, &_TMP21a);

        /* Measurement error e and covariance diag w */
        axis3f_t w;
        if (state.restDetected && params.restBiasEstEnabled) {
            _e.data[0] = state.restLastGyrLp.x - state.bias.x;
            _e.data[1] = state.restLastGyrLp.y - state.bias.y;
            _e.data[2] = state.restLastGyrLp.z - state.bias.z;
            vqf_setMatIdentityScaled(&_R, 1.0f);
            w.x = coeffs.biasRestW;
            w.y = coeffs.biasRestW;
            w.z = coeffs.biasRestW;
        } else if (params.motionBiasEstEnabled) {
            /* Recompute R*bias after filtering */
            matrixMult(&_R, &biasVec, &_TMP31a);
            _e.data[0] = -accEarth.y / coeffs.accTs + _TMP21a.data[0] - _TMP31a.data[0];
            _e.data[1] = accEarth.x / coeffs.accTs + _TMP21a.data[1] - _TMP31a.data[1];
            _e.data[2] = -_TMP31a.data[2];
            w.x = coeffs.biasMotionW;
            w.y = coeffs.biasMotionW;
            w.z = coeffs.biasVerticalW;
        } else {
            w.x = -1.0f;
            w.y = -1.0f;
            w.z = -1.0f;
        }

        /* Step 1: increase covariance even without update */
        if (ELEM(state.biasP, 0, 0) < coeffs.biasP0) {
            ELEM(state.biasP, 0, 0) += coeffs.biasV;
        }
        if (ELEM(state.biasP, 1, 1) < coeffs.biasP0) {
            ELEM(state.biasP, 1, 1) += coeffs.biasV;
        }
        if (ELEM(state.biasP, 2, 2) < coeffs.biasP0) {
            ELEM(state.biasP, 2, 2) += coeffs.biasV;
        }

        if (w.x >= 0.0f) {
            /* Clip disagreement */
            _e.data[0] = CONSTRAIN(_e.data[0], -biasClipRad, biasClipRad);
            _e.data[1] = CONSTRAIN(_e.data[1], -biasClipRad, biasClipRad);
            _e.data[2] = CONSTRAIN(_e.data[2], -biasClipRad, biasClipRad);

            /* Step 2: K = P R^T inv(W + R P R^T) */
            matrixMult_rhsT(&state.biasP, &_R, &_TMP33a); /* TMP33a = P R^T */
            matrixMult(&_R, &_TMP33a, &_TMP33b);          /* TMP33b = R P R^T */
            ELEM(_TMP33b, 0, 0) += w.x;
            ELEM(_TMP33b, 1, 1) += w.y;
            ELEM(_TMP33b, 2, 2) += w.z;
            matrixInversed(&_TMP33b, &_TMP33c);       /* TMP33c = inv(...) */
            matrixMult(&_TMP33a, &_TMP33c, &_TMP33d); /* TMP33d = K */

            /* Step 3: bias = bias + K e */
            matrixMult(&_TMP33d, &_e, &_TMP31b);
            state.bias.x += _TMP31b.data[0];
            state.bias.y += _TMP31b.data[1];
            state.bias.z += _TMP31b.data[2];

            /* Step 4: P = P - K R P */
            matrixMult(&_TMP33d, &_R, &_TMP33a);          /* TMP33a = K R */
            matrixMult(&_TMP33a, &state.biasP, &_TMP33b); /* TMP33b = K R P */
            matrixSub(&state.biasP, &_TMP33b, &state.biasP);

            /* Clip bias estimate */
            state.bias.x = CONSTRAIN(state.bias.x, -biasClipRad, biasClipRad);
            state.bias.y = CONSTRAIN(state.bias.y, -biasClipRad, biasClipRad);
            state.bias.z = CONSTRAIN(state.bias.z, -biasClipRad, biasClipRad);
        }
    }
}

void AHRS_VQF_UpdateMag(axis3f_t mag) {

    /* Convert from body-NED to internal body-ENU */
    const axis3f_t magBody = vqf_bodyNedToEnu(mag);

    if (magBody.x == 0.0f && magBody.y == 0.0f && magBody.z == 0.0f) {
        return;
    }

    quaternion_t accGyrQuat;
    vqf_getQuat6D_enu(&accGyrQuat);

    axis3f_t magEarth = vqf_quatRotateForward(&accGyrQuat, magBody);

    if (params.magDistRejectionEnabled) {
        const float norm = vqf_vecNorm(magEarth);
        state.magNormDip.data[0] = norm;
        state.magNormDip.data[1] = -asinf(CONSTRAIN(magEarth.z / fmaxf(norm, VQF_EPS), -1.0f, 1.0f));

        if (params.magCurrentTau > 0.0f) {
            vqf_filterVec(&state.magNormDip, params.magCurrentTau, coeffs.magTs, &coeffs.magNormDipLpB, &coeffs.magNormDipLpA, &state.magNormDipLpState,
                          &state.magNormDip);
        }

        /* Disturbance detection */
        const float refNorm = state.magRefNorm;
        const float refDip = state.magRefDip;

        if (refNorm > 0.0f && fabsf(state.magNormDip.data[0] - refNorm) < params.magNormTh * refNorm
            && fabsf(state.magNormDip.data[1] - refDip) < params.magDipTh * (constPI / 180.0f)) {
            state.magUndisturbedT += coeffs.magTs;
            if (state.magUndisturbedT >= params.magMinUndisturbedTime) {
                state.magDistDetected = 0u;
                state.magRefNorm += coeffs.kMagRef * (state.magNormDip.data[0] - state.magRefNorm);
                state.magRefDip += coeffs.kMagRef * (state.magNormDip.data[1] - state.magRefDip);
            }
        } else {
            state.magUndisturbedT = 0.0f;
            state.magDistDetected = 1u;
        }

        /* New field acceptance */
        if (state.magCandidateNorm > 0.0f && fabsf(state.magNormDip.data[0] - state.magCandidateNorm) < params.magNormTh * state.magCandidateNorm
            && fabsf(state.magNormDip.data[1] - state.magCandidateDip) < params.magDipTh * (constPI / 180.0f)) {
            if (vqf_vecNorm(state.restLastGyrLp) >= params.magNewMinGyr * (constPI / 180.0f)) {
                state.magCandidateT += coeffs.magTs;
            }

            state.magCandidateNorm += coeffs.kMagRef * (state.magNormDip.data[0] - state.magCandidateNorm);
            state.magCandidateDip += coeffs.kMagRef * (state.magNormDip.data[1] - state.magCandidateDip);

            if (state.magDistDetected
                && (state.magCandidateT >= params.magNewTime || (state.magRefNorm == 0.0f && state.magCandidateT >= params.magNewFirstTime))) {
                state.magRefNorm = state.magCandidateNorm;
                state.magRefDip = state.magCandidateDip;
                state.magDistDetected = 0u;
                state.magUndisturbedT = params.magMinUndisturbedTime;
            }
        } else {
            state.magCandidateT = 0.0f;
            state.magCandidateNorm = state.magNormDip.data[0];
            state.magCandidateDip = state.magNormDip.data[1];
        }
    }

    /* Disagreement angle based on current mag measurement (internal ENU) */
    state.lastMagDisAngle = atan2f(magEarth.x, magEarth.y) - state.delta;
    state.lastMagDisAngle = vqf_wrapPi(state.lastMagDisAngle);

    float k = coeffs.kMag;

    if (params.magDistRejectionEnabled) {
        if (state.magDistDetected) {
            if (state.magRejectT <= params.magMaxRejectionTime) {
                state.magRejectT += coeffs.magTs;
                k = 0.0f;
            } else {
                k /= params.magRejectionFactor;
            }
        } else {
            state.magRejectT = fmaxf(state.magRejectT - params.magRejectionFactor * coeffs.magTs, 0.0f);
        }
    }

    /* Ensure fast initial convergence */
    if (state.kMagInit != 0.0f) {
        if (k < state.kMagInit) {
            k = state.kMagInit;
        }

        state.kMagInit = state.kMagInit / (state.kMagInit + 1.0f);
        if (state.kMagInit * params.tauMag < coeffs.magTs) {
            state.kMagInit = 0.0f;
        }
    }

    state.delta += k * state.lastMagDisAngle;
    state.lastMagCorrAngularRate = k * state.lastMagDisAngle / coeffs.magTs;
    state.delta = vqf_wrapPi(state.delta);
}

/*----------------------------Get values from state vector----------------------------------*/

void AHRS_VQF_GetQuat3D(quaternion_t* q_out) {
    if (!q_out) {
        return;
    }
    quaternion_t q_enu;
    vqf_getQuat3D_enu(&q_enu);
    vqf_quatEnuToNed(&q_enu, q_out);
}

void AHRS_VQF_GetQuat6D(quaternion_t* q_out) {
    if (!q_out) {
        return;
    }
    quaternion_t q_enu;
    vqf_getQuat6D_enu(&q_enu);
    vqf_quatEnuToNed(&q_enu, q_out);
}

void AHRS_VQF_GetQuat9D(quaternion_t* q_out) {
    if (!q_out) {
        return;
    }
    quaternion_t q_enu;
    vqf_getQuat9D_enu(&q_enu);
    vqf_quatEnuToNed(&q_enu, q_out);
}

float AHRS_VQF_GetDelta() {
    /* Convert +Up (ENU) to +Down (NED) */
    return -state.delta;
}

float AHRS_VQF_GetBiasEstimate(axis3f_t* bias_out) {
    if (bias_out) {
        *bias_out = vqf_bodyEnuToNed(state.bias);
    }

    /* Upper bound on max eigenvalue (Gershgorin) and clip to biasP0 */
    const float sum1 = fabsf(ELEM(state.biasP, 0, 0)) + fabsf(ELEM(state.biasP, 0, 1)) + fabsf(ELEM(state.biasP, 0, 2));
    const float sum2 = fabsf(ELEM(state.biasP, 1, 0)) + fabsf(ELEM(state.biasP, 1, 1)) + fabsf(ELEM(state.biasP, 1, 2));
    const float sum3 = fabsf(ELEM(state.biasP, 2, 0)) + fabsf(ELEM(state.biasP, 2, 1)) + fabsf(ELEM(state.biasP, 2, 2));
    const float P = fminf(fmaxf(fmaxf(sum1, sum2), sum3), coeffs.biasP0);

    /* Convert std from 0.01 deg/s to rad/s */
    return SQRT(P) * (constPI / (100.0f * 180.0f));
}

void AHRS_VQF_SetBiasEstimate(axis3f_t bias, float sigma) {
    state.bias = vqf_bodyNedToEnu(bias);
    if (sigma > 0.0f) {
        const float P = vqf_square(sigma * (180.0f * 100.0f / constPI));
        vqf_setMatIdentityScaled(&state.biasP, P);
    }
}

uint8_t AHRS_VQF_GetRestDetected() { return state.restDetected; }

uint8_t AHRS_VQF_GetMagDistDetected() { return state.magDistDetected; }

void AHRS_VQF_SetMagRef(float norm, float dip_rad) {
    state.magRefNorm = norm;
    /* Internally dip is stored as positive down (see original VQF: -asin(z/norm)) */
    state.magRefDip = dip_rad;
}
