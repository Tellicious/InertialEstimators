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

/* Functions -----------------------------------------------------------------*/

static inline float vqf_square(float x) {
    return x * x;
}

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

    outA_2x1->data[0] = 2.0f * (C * C - 1.0f) / D;         /* a1 (a0 = 1) */
    outA_2x1->data[1] = (1.0f - sqrt2 * C + C * C) / D;    /* a2 */
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

static float vqf_vecNorm(axis3f_t v) {
    return SQRT(v.x * v.x + v.y * v.y + v.z * v.z);
}

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

static void vqf_getQuat3D_enu(const AHRS_VQF_t* vqf, quaternion_t* q_out) {
    *q_out = vqf->state.gyrQuat;
}

static void vqf_getQuat6D_enu(const AHRS_VQF_t* vqf, quaternion_t* q_out) {
    quaternionMult((quaternion_t*)&vqf->state.accQuat, (quaternion_t*)&vqf->state.gyrQuat, q_out);
}

static void vqf_getQuat9D_enu(const AHRS_VQF_t* vqf, quaternion_t* q_out) {
    vqf_getQuat6D_enu(vqf, q_out);
    vqf_applyDelta(q_out, vqf->state.delta);
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

static void vqf_updateCoeffs(AHRS_VQF_t* vqf) {
    AHRS_VQF_Params_t* p = &vqf->params;
    AHRS_VQF_Coeffs_t* c = &vqf->coeffs;

    vqf_filterCoeffs(p->tauAcc, c->accTs, &c->accLpB, &c->accLpA);
    c->kMag = vqf_gainFromTau(p->tauMag, c->magTs);

    /* Bias covariance is kept in units of (0.01 deg/s)^2 (same as original VQF) */
    c->biasP0 = vqf_square(p->biasSigmaInit * 100.0f);
    c->biasV = vqf_square(0.1f * 100.0f) * c->accTs / p->biasForgettingTime;

    const float pMotion = vqf_square(p->biasSigmaMotion * 100.0f);
    c->biasMotionW = vqf_square(pMotion) / c->biasV + pMotion;
    c->biasVerticalW = c->biasMotionW / MAX(p->biasVerticalForgettingFactor, 1e-10f);

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

static void vqf_resetState(AHRS_VQF_t* vqf) {
    /* Quaternions */
    vqf->state.gyrQuat.q0 = 1.0f;
    vqf->state.gyrQuat.q1 = 0.0f;
    vqf->state.gyrQuat.q2 = 0.0f;
    vqf->state.gyrQuat.q3 = 0.0f;

    vqf->state.accQuat.q0 = 1.0f;
    vqf->state.accQuat.q1 = 0.0f;
    vqf->state.accQuat.q2 = 0.0f;
    vqf->state.accQuat.q3 = 0.0f;

    vqf->state.delta = 0.0f;

    vqf->state.restDetected = 0u;
    vqf->state.magDistDetected = 1u;

    vqf->state.lastAccLp.x = 0.0f;
    vqf->state.lastAccLp.y = 0.0f;
    vqf->state.lastAccLp.z = 0.0f;
    vqf_matrixFill(&vqf->state.accLpState, NAN);
    vqf->state.lastAccCorrAngularRate = 0.0f;

    vqf->state.kMagInit = 1.0f;
    vqf->state.lastMagDisAngle = 0.0f;
    vqf->state.lastMagCorrAngularRate = 0.0f;

    vqf->state.bias.x = 0.0f;
    vqf->state.bias.y = 0.0f;
    vqf->state.bias.z = 0.0f;
    vqf_setMatIdentityScaled(&vqf->state.biasP, vqf->coeffs.biasP0);

    /* Motion bias estimator LPFs */
    vqf_matrixFill(&vqf->state.motionBiasEstRLpState, NAN);
    vqf_matrixFill(&vqf->state.motionBiasEstBiasLpState, NAN);

    /* Rest detection */
    vqf_matrixFill(&vqf->state.restLastSquaredDeviations, 0.0f);
    vqf->state.restT = 0.0f;
    vqf->state.restLastGyrLp.x = 0.0f;
    vqf->state.restLastGyrLp.y = 0.0f;
    vqf->state.restLastGyrLp.z = 0.0f;
    vqf_matrixFill(&vqf->state.restGyrLpState, NAN);
    vqf->state.restLastAccLp.x = 0.0f;
    vqf->state.restLastAccLp.y = 0.0f;
    vqf->state.restLastAccLp.z = 0.0f;
    vqf_matrixFill(&vqf->state.restAccLpState, NAN);

    /* Magnetic reference and disturbance rejection */
    vqf->state.magRefNorm = 0.0f;
    vqf->state.magRefDip = 0.0f;
    vqf->state.magUndisturbedT = 0.0f;
    vqf->state.magRejectT = vqf->params.magMaxRejectionTime;

    vqf->state.magCandidateNorm = -1.0f;
    vqf->state.magCandidateDip = 0.0f;
    vqf->state.magCandidateT = 0.0f;

    vqf_matrixFill(&vqf->state.magNormDip, 0.0f);
    vqf_matrixFill(&vqf->state.magNormDipLpState, NAN);
}

/*------------------------------------Public functions--------------------------------------------*/

/*------------------------------------Initialization--------------------------------------------*/

void AHRS_VQF_Init(AHRS_VQF_t* vqf, float gyrTs_s, float accTs_s, float magTs_s) {
    if (!vqf) {
        return;
    }

    memset(vqf, 0, sizeof(*vqf));
    vqf_initDefaultParams(&vqf->params);

    /* Sampling times */
    vqf->coeffs.gyrTs = gyrTs_s;
    vqf->coeffs.accTs = (accTs_s > 0.0f) ? accTs_s : gyrTs_s;
    vqf->coeffs.magTs = (magTs_s > 0.0f) ? magTs_s : gyrTs_s;

    /* Allocate coefficient matrices */
    matrixInit(&vqf->coeffs.accLpB, 3, 1);
    matrixInit(&vqf->coeffs.accLpA, 2, 1);
    matrixInit(&vqf->coeffs.restGyrLpB, 3, 1);
    matrixInit(&vqf->coeffs.restGyrLpA, 2, 1);
    matrixInit(&vqf->coeffs.restAccLpB, 3, 1);
    matrixInit(&vqf->coeffs.restAccLpA, 2, 1);
    matrixInit(&vqf->coeffs.magNormDipLpB, 3, 1);
    matrixInit(&vqf->coeffs.magNormDipLpA, 2, 1);

    /* Allocate state matrices */
    matrixInit(&vqf->state.accLpState, 3, 2);
    matrixInit(&vqf->state.biasP, 3, 3);
    matrixInit(&vqf->state.motionBiasEstRLpState, 9, 2);
    matrixInit(&vqf->state.motionBiasEstBiasLpState, 2, 2);
    matrixInit(&vqf->state.restLastSquaredDeviations, 2, 1);
    matrixInit(&vqf->state.restGyrLpState, 3, 2);
    matrixInit(&vqf->state.restAccLpState, 3, 2);
    matrixInit(&vqf->state.magNormDip, 2, 1);
    matrixInit(&vqf->state.magNormDipLpState, 2, 2);

    /* Allocate scratch matrices */
    matrixInit(&vqf->_R, 3, 3);
    matrixInit(&vqf->_TMP33a, 3, 3);
    matrixInit(&vqf->_TMP33b, 3, 3);
    matrixInit(&vqf->_TMP33c, 3, 3);
    matrixInit(&vqf->_TMP33d, 3, 3);
    matrixInit(&vqf->_TMP31a, 3, 1);
    matrixInit(&vqf->_TMP31b, 3, 1);
    matrixInit(&vqf->_TMP21a, 2, 1);
    matrixInit(&vqf->_e, 3, 1);

    vqf_updateCoeffs(vqf);
    vqf_resetState(vqf);
}

void AHRS_VQF_Deinit(AHRS_VQF_t* vqf) {
    if (!vqf) {
        return;
    }

    /* Coeff matrices */
    matrixDelete(&vqf->coeffs.accLpB);
    matrixDelete(&vqf->coeffs.accLpA);
    matrixDelete(&vqf->coeffs.restGyrLpB);
    matrixDelete(&vqf->coeffs.restGyrLpA);
    matrixDelete(&vqf->coeffs.restAccLpB);
    matrixDelete(&vqf->coeffs.restAccLpA);
    matrixDelete(&vqf->coeffs.magNormDipLpB);
    matrixDelete(&vqf->coeffs.magNormDipLpA);

    /* State matrices */
    matrixDelete(&vqf->state.accLpState);
    matrixDelete(&vqf->state.biasP);
    matrixDelete(&vqf->state.motionBiasEstRLpState);
    matrixDelete(&vqf->state.motionBiasEstBiasLpState);
    matrixDelete(&vqf->state.restLastSquaredDeviations);
    matrixDelete(&vqf->state.restGyrLpState);
    matrixDelete(&vqf->state.restAccLpState);
    matrixDelete(&vqf->state.magNormDip);
    matrixDelete(&vqf->state.magNormDipLpState);

    /* Scratch matrices */
    matrixDelete(&vqf->_R);
    matrixDelete(&vqf->_TMP33a);
    matrixDelete(&vqf->_TMP33b);
    matrixDelete(&vqf->_TMP33c);
    matrixDelete(&vqf->_TMP33d);
    matrixDelete(&vqf->_TMP31a);
    matrixDelete(&vqf->_TMP31b);
    matrixDelete(&vqf->_TMP21a);
    matrixDelete(&vqf->_e);

    memset(vqf, 0, sizeof(*vqf));
}

void AHRS_VQF_Reset(AHRS_VQF_t* vqf) {
    if (!vqf) {
        return;
    }
    vqf_updateCoeffs(vqf);
    vqf_resetState(vqf);
}


/*------------------------------------Update functions--------------------------------------------*/

void AHRS_VQF_UpdateGyr(AHRS_VQF_t* vqf, axis3f_t gyr_rad_s) {
    if (!vqf) {
        return;
    }

    /* Convert from body-NED to internal body-ENU */
    const axis3f_t gyr = vqf_bodyNedToEnu(gyr_rad_s);

    /* Rest detection based on gyro */
    if (vqf->params.restBiasEstEnabled || vqf->params.magDistRejectionEnabled) {
        matrix_t x = { (float*)&gyr, 3, 1 };
        matrix_t out = { (float*)&vqf->state.restLastGyrLp, 3, 1 };
        vqf_filterVec(&x, vqf->params.restFilterTau, vqf->coeffs.gyrTs, &vqf->coeffs.restGyrLpB, &vqf->coeffs.restGyrLpA, &vqf->state.restGyrLpState, &out);

        const float dx = gyr.x - vqf->state.restLastGyrLp.x;
        const float dy = gyr.y - vqf->state.restLastGyrLp.y;
        const float dz = gyr.z - vqf->state.restLastGyrLp.z;
        vqf->state.restLastSquaredDeviations.data[0] = dx * dx + dy * dy + dz * dz;

        const float biasClipRad = vqf->params.biasClip * (constPI / 180.0f);
        const float thGyrRad = vqf->params.restThGyr * (constPI / 180.0f);
        if (vqf->state.restLastSquaredDeviations.data[0] >= thGyrRad * thGyrRad ||
            fabsf(vqf->state.restLastGyrLp.x) > biasClipRad ||
            fabsf(vqf->state.restLastGyrLp.y) > biasClipRad ||
            fabsf(vqf->state.restLastGyrLp.z) > biasClipRad) {
            vqf->state.restT = 0.0f;
            vqf->state.restDetected = 0u;
        }
    }

    /* Remove estimated bias (internal ENU) */
    axis3f_t gyrNoBias;
    gyrNoBias.x = gyr.x - vqf->state.bias.x;
    gyrNoBias.y = gyr.y - vqf->state.bias.y;
    gyrNoBias.z = gyr.z - vqf->state.bias.z;

    /* Gyro integration */
    const float gyrNorm = vqf_vecNorm(gyrNoBias);
    const float angle = gyrNorm * vqf->coeffs.gyrTs;
    if (gyrNorm > VQF_EPS) {
        quaternion_t dq;
        const float half = 0.5f * angle;
        const float s = SIN(half) / gyrNorm;
        dq.q0 = COS(half);
        dq.q1 = s * gyrNoBias.x;
        dq.q2 = s * gyrNoBias.y;
        dq.q3 = s * gyrNoBias.z;
        quaternionMult(&vqf->state.gyrQuat, &dq, &vqf->state.gyrQuat);
        quaternionNorm(&vqf->state.gyrQuat);
    }
}

void AHRS_VQF_UpdateAcc(AHRS_VQF_t* vqf, axis3f_t acc) {
    if (!vqf) {
        return;
    }

    /* Convert from body-NED to internal body-ENU */
    const axis3f_t accBody = vqf_bodyNedToEnu(acc);

    /* Ignore [0 0 0] */
    if (accBody.x == 0.0f && accBody.y == 0.0f && accBody.z == 0.0f) {
        return;
    }

    /* Rest detection based on accelerometer */
    if (vqf->params.restBiasEstEnabled) {
        matrix_t x = { (float*)&accBody, 3, 1 };
        matrix_t out = { (float*)&vqf->state.restLastAccLp, 3, 1 };
        vqf_filterVec(&x, vqf->params.restFilterTau, vqf->coeffs.accTs, &vqf->coeffs.restAccLpB, &vqf->coeffs.restAccLpA, &vqf->state.restAccLpState, &out);

        const float dx = accBody.x - vqf->state.restLastAccLp.x;
        const float dy = accBody.y - vqf->state.restLastAccLp.y;
        const float dz = accBody.z - vqf->state.restLastAccLp.z;
        vqf->state.restLastSquaredDeviations.data[1] = dx * dx + dy * dy + dz * dz;

        if (vqf->state.restLastSquaredDeviations.data[1] >= vqf->params.restThAcc * vqf->params.restThAcc) {
            vqf->state.restT = 0.0f;
            vqf->state.restDetected = 0u;
        } else {
            vqf->state.restT += vqf->coeffs.accTs;
            if (vqf->state.restT >= vqf->params.restMinT) {
                vqf->state.restDetected = 1u;
            }
        }
    }

    /* Filter acc in inertial frame (internal ENU) */
    axis3f_t accEarth = vqf_quatRotateForward(&vqf->state.gyrQuat, accBody);

    matrix_t accEarthMat = { (float*)&accEarth, 3, 1 };
    matrix_t lastAccLpMat = { (float*)&vqf->state.lastAccLp, 3, 1 };
    vqf_filterVec(&accEarthMat, vqf->params.tauAcc, vqf->coeffs.accTs, &vqf->coeffs.accLpB, &vqf->coeffs.accLpA, &vqf->state.accLpState, &lastAccLpMat);

    /* Transform to 6D earth frame and normalize */
    accEarth = vqf_quatRotateForward(&vqf->state.accQuat, vqf->state.lastAccLp);
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

    quaternionMult(&accCorr, &vqf->state.accQuat, &vqf->state.accQuat);
    quaternionNorm(&vqf->state.accQuat);

    vqf->state.lastAccCorrAngularRate = acosf(CONSTRAIN(accEarth.z, -1.0f, 1.0f)) / vqf->coeffs.accTs;

    /* Gyro bias estimation */
    if (vqf->params.motionBiasEstEnabled || vqf->params.restBiasEstEnabled) {
        const float biasClipRad = vqf->params.biasClip * (constPI / 180.0f);

        quaternion_t accGyrQuat;
        vqf_getQuat6D_enu(vqf, &accGyrQuat);

        /* Build rotation matrix R(accGyrQuat) in internal ENU */
        const float q0 = accGyrQuat.q0;
        const float q1 = accGyrQuat.q1;
        const float q2 = accGyrQuat.q2;
        const float q3 = accGyrQuat.q3;

        ELEM(vqf->_R, 0, 0) = 1.0f - 2.0f * (q2 * q2 + q3 * q3);
        ELEM(vqf->_R, 0, 1) = 2.0f * (q2 * q1 - q0 * q3);
        ELEM(vqf->_R, 0, 2) = 2.0f * (q0 * q2 + q3 * q1);
        ELEM(vqf->_R, 1, 0) = 2.0f * (q0 * q3 + q2 * q1);
        ELEM(vqf->_R, 1, 1) = 1.0f - 2.0f * (q1 * q1 + q3 * q3);
        ELEM(vqf->_R, 1, 2) = 2.0f * (q2 * q3 - q1 * q0);
        ELEM(vqf->_R, 2, 0) = 2.0f * (q3 * q1 - q0 * q2);
        ELEM(vqf->_R, 2, 1) = 2.0f * (q0 * q1 + q3 * q2);
        ELEM(vqf->_R, 2, 2) = 1.0f - 2.0f * (q1 * q1 + q2 * q2);

        /* biasLp = (R * bias)_{x,y} */
        matrix_t biasVec = { (float*)&vqf->state.bias, 3, 1 };
        matrixMult(&vqf->_R, &biasVec, &vqf->_TMP31a);
        vqf->_TMP21a.data[0] = vqf->_TMP31a.data[0];
        vqf->_TMP21a.data[1] = vqf->_TMP31a.data[1];

        /* Low-pass filter R and biasLp */
        matrix_t Rvec = { vqf->_R.data, 9, 1 };
        vqf_filterVec(&Rvec, vqf->params.tauAcc, vqf->coeffs.accTs, &vqf->coeffs.accLpB, &vqf->coeffs.accLpA, &vqf->state.motionBiasEstRLpState, &Rvec);
        vqf_filterVec(&vqf->_TMP21a, vqf->params.tauAcc, vqf->coeffs.accTs, &vqf->coeffs.accLpB, &vqf->coeffs.accLpA, &vqf->state.motionBiasEstBiasLpState, &vqf->_TMP21a);

        /* Measurement error e and covariance diag w */
        axis3f_t w;
        if (vqf->state.restDetected && vqf->params.restBiasEstEnabled) {
            vqf->_e.data[0] = vqf->state.restLastGyrLp.x - vqf->state.bias.x;
            vqf->_e.data[1] = vqf->state.restLastGyrLp.y - vqf->state.bias.y;
            vqf->_e.data[2] = vqf->state.restLastGyrLp.z - vqf->state.bias.z;
            vqf_setMatIdentityScaled(&vqf->_R, 1.0f);
            w.x = vqf->coeffs.biasRestW;
            w.y = vqf->coeffs.biasRestW;
            w.z = vqf->coeffs.biasRestW;
        } else if (vqf->params.motionBiasEstEnabled) {
            /* Recompute R*bias after filtering */
            matrixMult(&vqf->_R, &biasVec, &vqf->_TMP31a);
            vqf->_e.data[0] = -accEarth.y / vqf->coeffs.accTs + vqf->_TMP21a.data[0] - vqf->_TMP31a.data[0];
            vqf->_e.data[1] = accEarth.x / vqf->coeffs.accTs + vqf->_TMP21a.data[1] - vqf->_TMP31a.data[1];
            vqf->_e.data[2] = -vqf->_TMP31a.data[2];
            w.x = vqf->coeffs.biasMotionW;
            w.y = vqf->coeffs.biasMotionW;
            w.z = vqf->coeffs.biasVerticalW;
        } else {
            w.x = -1.0f;
            w.y = -1.0f;
            w.z = -1.0f;
        }

        /* Step 1: increase covariance even without update */
        if (ELEM(vqf->state.biasP, 0, 0) < vqf->coeffs.biasP0) {
            ELEM(vqf->state.biasP, 0, 0) += vqf->coeffs.biasV;
        }
        if (ELEM(vqf->state.biasP, 1, 1) < vqf->coeffs.biasP0) {
            ELEM(vqf->state.biasP, 1, 1) += vqf->coeffs.biasV;
        }
        if (ELEM(vqf->state.biasP, 2, 2) < vqf->coeffs.biasP0) {
            ELEM(vqf->state.biasP, 2, 2) += vqf->coeffs.biasV;
        }

        if (w.x >= 0.0f) {
            /* Clip disagreement */
            vqf->_e.data[0] = CONSTRAIN(vqf->_e.data[0], -biasClipRad, biasClipRad);
            vqf->_e.data[1] = CONSTRAIN(vqf->_e.data[1], -biasClipRad, biasClipRad);
            vqf->_e.data[2] = CONSTRAIN(vqf->_e.data[2], -biasClipRad, biasClipRad);

            /* Step 2: K = P R^T inv(W + R P R^T) */
            matrixMult_rhsT(&vqf->state.biasP, &vqf->_R, &vqf->_TMP33a); /* TMP33a = P R^T */
            matrixMult(&vqf->_R, &vqf->_TMP33a, &vqf->_TMP33b);          /* TMP33b = R P R^T */
            ELEM(vqf->_TMP33b, 0, 0) += w.x;
            ELEM(vqf->_TMP33b, 1, 1) += w.y;
            ELEM(vqf->_TMP33b, 2, 2) += w.z;
            matrixInversed(&vqf->_TMP33b, &vqf->_TMP33c);               /* TMP33c = inv(...) */
            matrixMult(&vqf->_TMP33a, &vqf->_TMP33c, &vqf->_TMP33d);    /* TMP33d = K */

            /* Step 3: bias = bias + K e */
            matrixMult(&vqf->_TMP33d, &vqf->_e, &vqf->_TMP31b);
            vqf->state.bias.x += vqf->_TMP31b.data[0];
            vqf->state.bias.y += vqf->_TMP31b.data[1];
            vqf->state.bias.z += vqf->_TMP31b.data[2];

            /* Step 4: P = P - K R P */
            matrixMult(&vqf->_TMP33d, &vqf->_R, &vqf->_TMP33a);          /* TMP33a = K R */
            matrixMult(&vqf->_TMP33a, &vqf->state.biasP, &vqf->_TMP33b); /* TMP33b = K R P */
            matrixSub(&vqf->state.biasP, &vqf->_TMP33b, &vqf->state.biasP);

            /* Clip bias estimate */
            vqf->state.bias.x = CONSTRAIN(vqf->state.bias.x, -biasClipRad, biasClipRad);
            vqf->state.bias.y = CONSTRAIN(vqf->state.bias.y, -biasClipRad, biasClipRad);
            vqf->state.bias.z = CONSTRAIN(vqf->state.bias.z, -biasClipRad, biasClipRad);
        }
    }
}

void AHRS_VQF_UpdateMag(AHRS_VQF_t* vqf, axis3f_t mag) {
    if (!vqf) {
        return;
    }

    /* Convert from body-NED to internal body-ENU */
    const axis3f_t magBody = vqf_bodyNedToEnu(mag);

    if (magBody.x == 0.0f && magBody.y == 0.0f && magBody.z == 0.0f) {
        return;
    }

    quaternion_t accGyrQuat;
    vqf_getQuat6D_enu(vqf, &accGyrQuat);

    axis3f_t magEarth = vqf_quatRotateForward(&accGyrQuat, magBody);

    if (vqf->params.magDistRejectionEnabled) {
        const float norm = vqf_vecNorm(magEarth);
        vqf->state.magNormDip.data[0] = norm;
        vqf->state.magNormDip.data[1] = -asinf(CONSTRAIN(magEarth.z / MAX(norm, VQF_EPS), -1.0f, 1.0f));

        if (vqf->params.magCurrentTau > 0.0f) {
            vqf_filterVec(&vqf->state.magNormDip, vqf->params.magCurrentTau, vqf->coeffs.magTs, &vqf->coeffs.magNormDipLpB, &vqf->coeffs.magNormDipLpA, &vqf->state.magNormDipLpState, &vqf->state.magNormDip);
        }

        /* Disturbance detection */
        const float refNorm = vqf->state.magRefNorm;
        const float refDip = vqf->state.magRefDip;

        if (refNorm > 0.0f &&
            fabsf(vqf->state.magNormDip.data[0] - refNorm) < vqf->params.magNormTh * refNorm &&
            fabsf(vqf->state.magNormDip.data[1] - refDip) < vqf->params.magDipTh * (constPI / 180.0f)) {
            vqf->state.magUndisturbedT += vqf->coeffs.magTs;
            if (vqf->state.magUndisturbedT >= vqf->params.magMinUndisturbedTime) {
                vqf->state.magDistDetected = 0u;
                vqf->state.magRefNorm += vqf->coeffs.kMagRef * (vqf->state.magNormDip.data[0] - vqf->state.magRefNorm);
                vqf->state.magRefDip += vqf->coeffs.kMagRef * (vqf->state.magNormDip.data[1] - vqf->state.magRefDip);
            }
        } else {
            vqf->state.magUndisturbedT = 0.0f;
            vqf->state.magDistDetected = 1u;
        }

        /* New field acceptance */
        if (vqf->state.magCandidateNorm > 0.0f &&
            fabsf(vqf->state.magNormDip.data[0] - vqf->state.magCandidateNorm) < vqf->params.magNormTh * vqf->state.magCandidateNorm &&
            fabsf(vqf->state.magNormDip.data[1] - vqf->state.magCandidateDip) < vqf->params.magDipTh * (constPI / 180.0f)) {
            if (vqf_vecNorm(vqf->state.restLastGyrLp) >= vqf->params.magNewMinGyr * (constPI / 180.0f)) {
                vqf->state.magCandidateT += vqf->coeffs.magTs;
            }

            vqf->state.magCandidateNorm += vqf->coeffs.kMagRef * (vqf->state.magNormDip.data[0] - vqf->state.magCandidateNorm);
            vqf->state.magCandidateDip += vqf->coeffs.kMagRef * (vqf->state.magNormDip.data[1] - vqf->state.magCandidateDip);

            if (vqf->state.magDistDetected &&
                (vqf->state.magCandidateT >= vqf->params.magNewTime || (vqf->state.magRefNorm == 0.0f && vqf->state.magCandidateT >= vqf->params.magNewFirstTime))) {
                vqf->state.magRefNorm = vqf->state.magCandidateNorm;
                vqf->state.magRefDip = vqf->state.magCandidateDip;
                vqf->state.magDistDetected = 0u;
                vqf->state.magUndisturbedT = vqf->params.magMinUndisturbedTime;
            }
        } else {
            vqf->state.magCandidateT = 0.0f;
            vqf->state.magCandidateNorm = vqf->state.magNormDip.data[0];
            vqf->state.magCandidateDip = vqf->state.magNormDip.data[1];
        }
    }

    /* Disagreement angle based on current mag measurement (internal ENU) */
    vqf->state.lastMagDisAngle = atan2f(magEarth.x, magEarth.y) - vqf->state.delta;
    vqf->state.lastMagDisAngle = vqf_wrapPi(vqf->state.lastMagDisAngle);

    float k = vqf->coeffs.kMag;

    if (vqf->params.magDistRejectionEnabled) {
        if (vqf->state.magDistDetected) {
            if (vqf->state.magRejectT <= vqf->params.magMaxRejectionTime) {
                vqf->state.magRejectT += vqf->coeffs.magTs;
                k = 0.0f;
            } else {
                k /= vqf->params.magRejectionFactor;
            }
        } else {
            vqf->state.magRejectT = MAX(vqf->state.magRejectT - vqf->params.magRejectionFactor * vqf->coeffs.magTs, 0.0f);
        }
    }

    /* Ensure fast initial convergence */
    if (vqf->state.kMagInit != 0.0f) {
        if (k < vqf->state.kMagInit) {
            k = vqf->state.kMagInit;
        }

        vqf->state.kMagInit = vqf->state.kMagInit / (vqf->state.kMagInit + 1.0f);
        if (vqf->state.kMagInit * vqf->params.tauMag < vqf->coeffs.magTs) {
            vqf->state.kMagInit = 0.0f;
        }
    }

    vqf->state.delta += k * vqf->state.lastMagDisAngle;
    vqf->state.lastMagCorrAngularRate = k * vqf->state.lastMagDisAngle / vqf->coeffs.magTs;
    vqf->state.delta = vqf_wrapPi(vqf->state.delta);
}


/*----------------------------Get values from state vector----------------------------------*/

void AHRS_VQF_GetQuat3D(const AHRS_VQF_t* vqf, quaternion_t* q_out) {
    if (!vqf || !q_out) {
        return;
    }
    quaternion_t q_enu;
    vqf_getQuat3D_enu(vqf, &q_enu);
    vqf_quatEnuToNed(&q_enu, q_out);
}

void AHRS_VQF_GetQuat6D(const AHRS_VQF_t* vqf, quaternion_t* q_out) {
    if (!vqf || !q_out) {
        return;
    }
    quaternion_t q_enu;
    vqf_getQuat6D_enu(vqf, &q_enu);
    vqf_quatEnuToNed(&q_enu, q_out);
}

void AHRS_VQF_GetQuat9D(const AHRS_VQF_t* vqf, quaternion_t* q_out) {
    if (!vqf || !q_out) {
        return;
    }
    quaternion_t q_enu;
    vqf_getQuat9D_enu(vqf, &q_enu);
    vqf_quatEnuToNed(&q_enu, q_out);
}

float AHRS_VQF_GetDelta(const AHRS_VQF_t* vqf) {
    if (!vqf) {
        return 0.0f;
    }
    /* Convert +Up (ENU) to +Down (NED) */
    return -vqf->state.delta;
}

float AHRS_VQF_GetBiasEstimate(const AHRS_VQF_t* vqf, axis3f_t* bias_out) {
    if (!vqf) {
        return 0.0f;
    }

    if (bias_out) {
        *bias_out = vqf_bodyEnuToNed(vqf->state.bias);
    }

    /* Upper bound on max eigenvalue (Gershgorin) and clip to biasP0 */
    const float sum1 = fabsf(ELEM(vqf->state.biasP, 0, 0)) + fabsf(ELEM(vqf->state.biasP, 0, 1)) + fabsf(ELEM(vqf->state.biasP, 0, 2));
    const float sum2 = fabsf(ELEM(vqf->state.biasP, 1, 0)) + fabsf(ELEM(vqf->state.biasP, 1, 1)) + fabsf(ELEM(vqf->state.biasP, 1, 2));
    const float sum3 = fabsf(ELEM(vqf->state.biasP, 2, 0)) + fabsf(ELEM(vqf->state.biasP, 2, 1)) + fabsf(ELEM(vqf->state.biasP, 2, 2));
    const float P = MIN(MAX(MAX(sum1, sum2), sum3), vqf->coeffs.biasP0);

    /* Convert std from 0.01 deg/s to rad/s */
    return SQRT(P) * (constPI / (100.0f * 180.0f));
}

void AHRS_VQF_SetBiasEstimate(AHRS_VQF_t* vqf, axis3f_t bias, float sigma) {
    if (!vqf) {
        return;
    }

    vqf->state.bias = vqf_bodyNedToEnu(bias);
    if (sigma > 0.0f) {
        const float P = vqf_square(sigma * (180.0f * 100.0f / constPI));
        vqf_setMatIdentityScaled(&vqf->state.biasP, P);
    }
}

uint8_t AHRS_VQF_GetRestDetected(const AHRS_VQF_t* vqf) {
    if (!vqf) {
        return 0u;
    }
    return vqf->state.restDetected;
}

uint8_t AHRS_VQF_GetMagDistDetected(const AHRS_VQF_t* vqf) {
    if (!vqf) {
        return 0u;
    }
    return vqf->state.magDistDetected;
}

void AHRS_VQF_SetMagRef(AHRS_VQF_t* vqf, float norm, float dip_rad) {
    if (!vqf) {
        return;
    }
    vqf->state.magRefNorm = norm;
    /* Internally dip is stored as positive down (see original VQF: -asin(z/norm)) */
    vqf->state.magRefDip = dip_rad;
}
