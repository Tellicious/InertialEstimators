/* BEGIN Header */
/**
 ******************************************************************************
 * \file            AHRS_VQF.c
 * \author          Andrea Vivani
 * \brief           VQF (Versatile Quaternion-based Filter) AHRS implementation
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

/* Includes ------------------------------------------------------------------*/

#include "AHRS_VQF.h"

#include <float.h>
#include <math.h>
#include <string.h>

/* Private defines -----------------------------------------------------------*/

#define VQF_EPS    (FLT_EPSILON)
#define VQF_NAN    (NAN)
#define VQF_SQRT2  (1.41421356237309504880f)

/* Private function prototypes -----------------------------------------------*/

static void vqf_bind_matrices(AHRS_VQF_t* vqf);

static inline axis3f_t vqf_body_ned_to_enu(axis3f_t v_ned);
static inline axis3f_t vqf_body_enu_to_ned(axis3f_t v_enu);

static void vqf_quat_enu_to_ned(const quaternion_t* q_enu, quaternion_t* q_ned);
static void vqf_quat_apply_delta(const quaternion_t* q_in, float delta, quaternion_t* q_out);
static void vqf_quat_rotate(const quaternion_t* q, const axis3f_t* v, axis3f_t* out);

static float vqf_norm_f32(const float* vec, uint16_t n);
static void  vqf_normalize_f32(float* vec, uint16_t n);

static float vqf_gain_from_tau(float tau, float Ts);
static void  vqf_filter_coeffs(float tau, float Ts, float b[3], float a[2]);
static void  vqf_filter_initial_state(float x0, const float b[3], const float a[2], float state[2]);
static void  vqf_filter_adapt_state_for_coeff_change(const float* last_y, uint16_t n,
                                                     const float b_old[3], const float a_old[2],
                                                     const float b_new[3], const float a_new[2],
                                                     float* state);
static float vqf_filter_step(float x, const float b[3], const float a[2], float state[2]);
static void  vqf_filter_vec(const float* x, uint16_t n, float tau, float Ts,
                            const float b[3], const float a[2],
                            float* state, float* out);

static void  vqf_init_params(AHRS_VQF_Params_t* p);
static void  vqf_reset_state(AHRS_VQF_t* vqf);
static void  vqf_setup(AHRS_VQF_t* vqf);

static void  vqf_get_quat3d_enu(const AHRS_VQF_t* vqf, quaternion_t* q_out);
static void  vqf_get_quat6d_enu(const AHRS_VQF_t* vqf, quaternion_t* q_out);
static void  vqf_get_quat9d_enu(const AHRS_VQF_t* vqf, quaternion_t* q_out);

static void  vqf_mat3_set_scaled_identity(float scale, matrix_t* M);
static void  vqf_mat3_inv(const matrix_t* A, matrix_t* Ainv);

/* Private functions ---------------------------------------------------------*/

static void vqf_bind_matrices(AHRS_VQF_t* vqf)
{
    /* Coefficients */
    vqf->coeffs.accLpB.rows = 3u;
    vqf->coeffs.accLpB.cols = 1u;
    vqf->coeffs.accLpB.data = vqf->coeffs.accLpBData;

    vqf->coeffs.accLpA.rows = 2u;
    vqf->coeffs.accLpA.cols = 1u;
    vqf->coeffs.accLpA.data = vqf->coeffs.accLpAData;

    vqf->coeffs.restGyrLpB.rows = 3u;
    vqf->coeffs.restGyrLpB.cols = 1u;
    vqf->coeffs.restGyrLpB.data = vqf->coeffs.restGyrLpBData;

    vqf->coeffs.restGyrLpA.rows = 2u;
    vqf->coeffs.restGyrLpA.cols = 1u;
    vqf->coeffs.restGyrLpA.data = vqf->coeffs.restGyrLpAData;

    vqf->coeffs.restAccLpB.rows = 3u;
    vqf->coeffs.restAccLpB.cols = 1u;
    vqf->coeffs.restAccLpB.data = vqf->coeffs.restAccLpBData;

    vqf->coeffs.restAccLpA.rows = 2u;
    vqf->coeffs.restAccLpA.cols = 1u;
    vqf->coeffs.restAccLpA.data = vqf->coeffs.restAccLpAData;

    vqf->coeffs.magNormDipLpB.rows = 3u;
    vqf->coeffs.magNormDipLpB.cols = 1u;
    vqf->coeffs.magNormDipLpB.data = vqf->coeffs.magNormDipLpBData;

    vqf->coeffs.magNormDipLpA.rows = 2u;
    vqf->coeffs.magNormDipLpA.cols = 1u;
    vqf->coeffs.magNormDipLpA.data = vqf->coeffs.magNormDipLpAData;

    /* State */
    vqf->state.accLpState.rows = 6u;
    vqf->state.accLpState.cols = 1u;
    vqf->state.accLpState.data = vqf->state.accLpStateData;

    vqf->state.biasP.rows = 3u;
    vqf->state.biasP.cols = 3u;
    vqf->state.biasP.data = vqf->state.biasPData;

    vqf->state.motionBiasEstRLpState.rows = 18u;
    vqf->state.motionBiasEstRLpState.cols = 1u;
    vqf->state.motionBiasEstRLpState.data = vqf->state.motionBiasEstRLpStateData;

    vqf->state.motionBiasEstBiasLpState.rows = 4u;
    vqf->state.motionBiasEstBiasLpState.cols = 1u;
    vqf->state.motionBiasEstBiasLpState.data = vqf->state.motionBiasEstBiasLpStateData;

    vqf->state.restLastSquaredDeviations.rows = 2u;
    vqf->state.restLastSquaredDeviations.cols = 1u;
    vqf->state.restLastSquaredDeviations.data = vqf->state.restLastSquaredDeviationsData;

    vqf->state.restGyrLpState.rows = 6u;
    vqf->state.restGyrLpState.cols = 1u;
    vqf->state.restGyrLpState.data = vqf->state.restGyrLpStateData;

    vqf->state.restAccLpState.rows = 6u;
    vqf->state.restAccLpState.cols = 1u;
    vqf->state.restAccLpState.data = vqf->state.restAccLpStateData;

    vqf->state.magNormDip.rows = 2u;
    vqf->state.magNormDip.cols = 1u;
    vqf->state.magNormDip.data = vqf->state.magNormDipData;

    vqf->state.magNormDipLpState.rows = 4u;
    vqf->state.magNormDipLpState.cols = 1u;
    vqf->state.magNormDipLpState.data = vqf->state.magNormDipLpStateData;
}

static inline axis3f_t vqf_body_ned_to_enu(axis3f_t v_ned)
{
    /* body-NED -> body-ENU: [E, N, U] = [Y, X, -Z] */
    axis3f_t v_enu;
    v_enu.x = v_ned.y;
    v_enu.y = v_ned.x;
    v_enu.z = -v_ned.z;
    return v_enu;
}

static inline axis3f_t vqf_body_enu_to_ned(axis3f_t v_enu)
{
    /* body-ENU -> body-NED: [N, E, D] = [Y, X, -Z] */
    axis3f_t v_ned;
    v_ned.x = v_enu.y;
    v_ned.y = v_enu.x;
    v_ned.z = -v_enu.z;
    return v_ned;
}

static void vqf_quat_enu_to_ned(const quaternion_t* q_enu, quaternion_t* q_ned)
{
    /* Convert both body and navigation frames: q_ned = qP ⊗ q_enu ⊗ qP
     * where qP maps NED -> ENU (swap x/y and flip z). qP is self-inverse up to sign.
     */
    const float s = 0.70710678118654752440f;
    const quaternion_t qP = {0.0f, s, s, 0.0f};

    quaternion_t tmp;
    quaternionMult((quaternion_t*)&qP, (quaternion_t*)q_enu, &tmp);
    quaternionMult(&tmp, (quaternion_t*)&qP, q_ned);
    quaternionNorm(q_ned);
}

static void vqf_quat_apply_delta(const quaternion_t* q_in, float delta, quaternion_t* q_out)
{
    /* q_out = q_delta ⊗ q_in, q_delta is a rotation about +Z of the ENU navigation frame. */
    quaternion_t q_delta;
    q_delta.q0 = COS(0.5f * delta);
    q_delta.q1 = 0.0f;
    q_delta.q2 = 0.0f;
    q_delta.q3 = SIN(0.5f * delta);

    quaternionMult(&q_delta, (quaternion_t*)q_in, q_out);
    quaternionNorm(q_out);
}

static void vqf_quat_rotate(const quaternion_t* q, const axis3f_t* v, axis3f_t* out)
{
    /* Active rotation: out = R(q) * v */
    const float w = q->q0;
    const float x = q->q1;
    const float y = q->q2;
    const float z = q->q3;

    const float vx = v->x;
    const float vy = v->y;
    const float vz = v->z;

    out->x = (1.0f - 2.0f * y * y - 2.0f * z * z) * vx
           + 2.0f * (x * y - w * z) * vy
           + 2.0f * (x * z + w * y) * vz;

    out->y = 2.0f * (x * y + w * z) * vx
           + (1.0f - 2.0f * x * x - 2.0f * z * z) * vy
           + 2.0f * (y * z - w * x) * vz;

    out->z = 2.0f * (x * z - w * y) * vx
           + 2.0f * (y * z + w * x) * vy
           + (1.0f - 2.0f * x * x - 2.0f * y * y) * vz;
}

static float vqf_norm_f32(const float* vec, uint16_t n)
{
    float s = 0.0f;
    for (uint16_t i = 0; i < n; i++) {
        s += vec[i] * vec[i];
    }
    return SQRT(s);
}

static void vqf_normalize_f32(float* vec, uint16_t n)
{
    const float nrm = vqf_norm_f32(vec, n);
    if (nrm <= VQF_EPS) {
        return;
    }
    const float inv = 1.0f / nrm;
    for (uint16_t i = 0; i < n; i++) {
        vec[i] *= inv;
    }
}

static float vqf_gain_from_tau(float tau, float Ts)
{
    if (tau < 0.0f) {
        return 0.0f;
    } else if (tau == 0.0f) {
        return 1.0f;
    } else {
        return 1.0f - expf(-Ts / tau);
    }
}

static void vqf_filter_coeffs(float tau, float Ts, float b[3], float a[2])
{
    /* 2nd order Butterworth filter based on https://stackoverflow.com/a/52764064
     * time constant of damped, non-oscillating part of step response.
     */
    if (tau <= 0.0f) {
        /* Tau==0 -> no filtering (pass-through) */
        b[0] = 1.0f; b[1] = 0.0f; b[2] = 0.0f;
        a[0] = 0.0f; a[1] = 0.0f;
        return;
    }

    const float fc = (VQF_SQRT2 / (2.0f * constPI)) / tau;
    const float C = TAN(constPI * fc * Ts);
    const float D = C * C + VQF_SQRT2 * C + 1.0f;

    const float b0 = (C * C) / D;
    b[0] = b0;
    b[1] = 2.0f * b0;
    b[2] = b0;

    a[0] = 2.0f * (C * C - 1.0f) / D;
    a[1] = (1.0f - VQF_SQRT2 * C + C * C) / D;
}

static void vqf_filter_initial_state(float x0, const float b[3], const float a[2], float state[2])
{
    /* initial state for steady-state (equilibrium) with constant input x0 */
    state[0] = x0 * (1.0f - b[0]);
    state[1] = (b[2] - a[1]) * x0;
}

static void vqf_filter_adapt_state_for_coeff_change(const float* last_y, uint16_t n,
                                                    const float b_old[3], const float a_old[2],
                                                    const float b_new[3], const float a_new[2],
                                                    float* state)
{
    /* Adapt filter state in-place to avoid discontinuity when coefficients change. */
    if (isnan(state[0])) {
        return; /* still in initialization phase */
    }

    for (uint16_t i = 0; i < n; i++) {
        state[2u * i + 0u] += (b_old[0] - b_new[0]) * last_y[i];
        state[2u * i + 1u] += (b_old[1] - b_new[1] - a_old[0] + a_new[0]) * last_y[i];
    }
}

static float vqf_filter_step(float x, const float b[3], const float a[2], float state[2])
{
    /* Direct Form II transposed */
    float y = state[0] + b[0] * x;

    state[0] = state[1] + b[1] * x - a[0] * y;
    state[1] = b[2] * x - a[1] * y;

    return y;
}

static void vqf_filter_vec(const float* x, uint16_t n, float tau, float Ts,
                           const float b[3], const float a[2],
                           float* state, float* out)
{
    /* state layout is identical to the reference implementation:
     * - During initialization phase: state[0] is NaN, state[1] holds sample count,
     *   and state[2+i] accumulates a running sum of the i-th signal.
     * - After initialization: state contains 2 states per signal (length n*2).
     */
    if (isnan(state[0])) {
        if (isnan(state[1])) {
            /* first sample */
            state[1] = 0.0f;
            for (uint16_t i = 0; i < n; i++) {
                state[2u + i] = 0.0f;
            }
        }

        state[1] += 1.0f;
        for (uint16_t i = 0; i < n; i++) {
            state[2u + i] += x[i];
            out[i] = state[2u + i] / state[1];
        }

        if (state[1] * Ts >= tau) {
            /* initialize the real filter state based on the averaged output */
            for (uint16_t i = 0; i < n; i++) {
                vqf_filter_initial_state(out[i], b, a, &state[2u * i]);
            }
        }
        return;
    }

    for (uint16_t i = 0; i < n; i++) {
        out[i] = vqf_filter_step(x[i], b, a, &state[2u * i]);
    }
}

/* Matrix helpers ------------------------------------------------------------*/

static void vqf_mat3_set_scaled_identity(float scale, matrix_t* M)
{
    if (M == NULL || M->data == NULL || M->rows != 3u || M->cols != 3u) {
        return;
    }

    matrixZeros(M);
    M->data[0] = scale;
    M->data[4] = scale;
    M->data[8] = scale;
}

static void vqf_mat3_inv(const matrix_t* A, matrix_t* Ainv)
{
    if (A == NULL || Ainv == NULL || A->data == NULL || Ainv->data == NULL ||
        A->rows != 3u || A->cols != 3u || Ainv->rows != 3u || Ainv->cols != 3u) {
        return;
    }

    const float a00 = A->data[0], a01 = A->data[1], a02 = A->data[2];
    const float a10 = A->data[3], a11 = A->data[4], a12 = A->data[5];
    const float a20 = A->data[6], a21 = A->data[7], a22 = A->data[8];

    const float det = a00 * (a11 * a22 - a12 * a21)
                    - a01 * (a10 * a22 - a12 * a20)
                    + a02 * (a10 * a21 - a11 * a20);

    if (ABS(det) <= VQF_EPS) {
        vqf_mat3_set_scaled_identity(1.0f, Ainv);
        return;
    }

    const float invDet = 1.0f / det;

    Ainv->data[0] =  (a11 * a22 - a12 * a21) * invDet;
    Ainv->data[1] = -(a01 * a22 - a02 * a21) * invDet;
    Ainv->data[2] =  (a01 * a12 - a02 * a11) * invDet;

    Ainv->data[3] = -(a10 * a22 - a12 * a20) * invDet;
    Ainv->data[4] =  (a00 * a22 - a02 * a20) * invDet;
    Ainv->data[5] = -(a00 * a12 - a02 * a10) * invDet;

    Ainv->data[6] =  (a10 * a21 - a11 * a20) * invDet;
    Ainv->data[7] = -(a00 * a21 - a01 * a20) * invDet;
    Ainv->data[8] =  (a00 * a11 - a01 * a10) * invDet;
}

/* Parameter / setup ---------------------------------------------------------*/

static void vqf_init_params(AHRS_VQF_Params_t* p)
{
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

static void vqf_reset_state(AHRS_VQF_t* vqf)
{
    AHRS_VQF_State_t* s = &vqf->state;

    /* Quaternion identities */
    s->gyrQuat.q0 = 1.0f; s->gyrQuat.q1 = 0.0f; s->gyrQuat.q2 = 0.0f; s->gyrQuat.q3 = 0.0f;
    s->accQuat.q0 = 1.0f; s->accQuat.q1 = 0.0f; s->accQuat.q2 = 0.0f; s->accQuat.q3 = 0.0f;
    s->delta = 0.0f;

    s->restDetected = 0u;
    s->magDistDetected = 1u;

    s->lastAccLp.x = 0.0f; s->lastAccLp.y = 0.0f; s->lastAccLp.z = 0.0f;

    for (uint16_t i = 0; i < 6u; i++) {
        s->accLpState.data[i] = VQF_NAN;
    }
    s->lastAccCorrAngularRate = 0.0f;

    s->kMagInit = 1.0f;
    s->lastMagDisAngle = 0.0f;
    s->lastMagCorrAngularRate = 0.0f;

    s->bias.x = 0.0f; s->bias.y = 0.0f; s->bias.z = 0.0f;
    vqf_mat3_set_scaled_identity(vqf->coeffs.biasP0, &s->biasP);

    for (uint16_t i = 0; i < 18u; i++) {
        s->motionBiasEstRLpState.data[i] = VQF_NAN;
    }
    for (uint16_t i = 0; i < 4u; i++) {
        s->motionBiasEstBiasLpState.data[i] = VQF_NAN;
    }

    s->restLastSquaredDeviations.data[0] = 0.0f;
    s->restLastSquaredDeviations.data[1] = 0.0f;
    s->restT = 0.0f;

    s->restLastGyrLp.x = 0.0f; s->restLastGyrLp.y = 0.0f; s->restLastGyrLp.z = 0.0f;
    for (uint16_t i = 0; i < 6u; i++) {
        s->restGyrLpState.data[i] = VQF_NAN;
    }

    s->restLastAccLp.x = 0.0f; s->restLastAccLp.y = 0.0f; s->restLastAccLp.z = 0.0f;
    for (uint16_t i = 0; i < 6u; i++) {
        s->restAccLpState.data[i] = VQF_NAN;
    }

    s->magRefNorm = 0.0f;
    s->magRefDip = 0.0f;
    s->magUndisturbedT = 0.0f;
    s->magRejectT = vqf->params.magMaxRejectionTime;

    s->magCandidateNorm = -1.0f;
    s->magCandidateDip = 0.0f;
    s->magCandidateT = 0.0f;

    s->magNormDip.data[0] = 0.0f;
    s->magNormDip.data[1] = 0.0f;
    for (uint16_t i = 0; i < 4u; i++) {
        s->magNormDipLpState.data[i] = VQF_NAN;
    }
}

static void vqf_setup(AHRS_VQF_t* vqf)
{
    AHRS_VQF_Params_t* p = &vqf->params;
    AHRS_VQF_Coeffs_t* c = &vqf->coeffs;

    vqf_filter_coeffs(p->tauAcc, c->accTs, c->accLpB.data, c->accLpA.data);
    c->kMag = vqf_gain_from_tau(p->tauMag, c->magTs);

    c->biasP0 = (p->biasSigmaInit * 100.0f) * (p->biasSigmaInit * 100.0f);
    c->biasV  = (0.1f * 100.0f) * (0.1f * 100.0f) * (c->accTs / p->biasForgettingTime);

    const float pMotion = (p->biasSigmaMotion * 100.0f) * (p->biasSigmaMotion * 100.0f);
    c->biasMotionW = (pMotion * pMotion) / c->biasV + pMotion;
    const float vff = MAX(p->biasVerticalForgettingFactor, 1e-10f);
    c->biasVerticalW = c->biasMotionW / vff;

    const float pRest = (p->biasSigmaRest * 100.0f) * (p->biasSigmaRest * 100.0f);
    c->biasRestW = (pRest * pRest) / c->biasV + pRest;

    vqf_filter_coeffs(p->restFilterTau, c->gyrTs, c->restGyrLpB.data, c->restGyrLpA.data);
    vqf_filter_coeffs(p->restFilterTau, c->accTs, c->restAccLpB.data, c->restAccLpA.data);

    c->kMagRef = vqf_gain_from_tau(p->magRefTau, c->magTs);

    if (p->magCurrentTau > 0.0f) {
        vqf_filter_coeffs(p->magCurrentTau, c->magTs, c->magNormDipLpB.data, c->magNormDipLpA.data);
    } else {
        /* Mark coefficients as unused (NaN) */
        c->magNormDipLpB.data[0] = VQF_NAN;
        c->magNormDipLpB.data[1] = VQF_NAN;
        c->magNormDipLpB.data[2] = VQF_NAN;
        c->magNormDipLpA.data[0] = VQF_NAN;
        c->magNormDipLpA.data[1] = VQF_NAN;
    }

    vqf_reset_state(vqf);
}

/* Internal quaternion getters (ENU) -----------------------------------------*/

static void vqf_get_quat3d_enu(const AHRS_VQF_t* vqf, quaternion_t* q_out)
{
    *q_out = vqf->state.gyrQuat;
}

static void vqf_get_quat6d_enu(const AHRS_VQF_t* vqf, quaternion_t* q_out)
{
    quaternionMult((quaternion_t*)&vqf->state.accQuat, (quaternion_t*)&vqf->state.gyrQuat, q_out);
    quaternionNorm(q_out);
}

static void vqf_get_quat9d_enu(const AHRS_VQF_t* vqf, quaternion_t* q_out)
{
    quaternion_t q6;
    vqf_get_quat6d_enu(vqf, &q6);
    vqf_quat_apply_delta(&q6, vqf->state.delta, q_out);
}

/* Public API ----------------------------------------------------------------*/

void AHRS_VQF_Init(AHRS_VQF_t* vqf, float gyrTs_s, float accTs_s, float magTs_s)
{
    if (vqf == NULL) {
        return;
    }

    vqf_init_params(&vqf->params);
    vqf_bind_matrices(vqf);

    vqf->coeffs.gyrTs = gyrTs_s;
    vqf->coeffs.accTs = (accTs_s > 0.0f) ? accTs_s : gyrTs_s;
    vqf->coeffs.magTs = (magTs_s > 0.0f) ? magTs_s : gyrTs_s;

    vqf_setup(vqf);
}

void AHRS_VQF_Reset(AHRS_VQF_t* vqf)
{
    if (vqf == NULL) {
        return;
    }
    vqf_reset_state(vqf);
}

void AHRS_VQF_UpdateGyr(AHRS_VQF_t* vqf, axis3f_t gyr_rad_s)
{
    if (vqf == NULL) {
        return;
    }

    AHRS_VQF_Params_t* p = &vqf->params;
    AHRS_VQF_Coeffs_t* c = &vqf->coeffs;
    AHRS_VQF_State_t*  s = &vqf->state;

    /* Convert body-NED sample to internal body-ENU */
    const axis3f_t gyr = vqf_body_ned_to_enu(gyr_rad_s);

    /* Rest detection (gyro) */
    if (p->restBiasEstEnabled || p->magDistRejectionEnabled) {
        const float gyrVec[3] = {gyr.x, gyr.y, gyr.z};
        float gyrLp[3];

        vqf_filter_vec(gyrVec, 3u, p->restFilterTau, c->gyrTs,
                       c->restGyrLpB.data, c->restGyrLpA.data,
                       s->restGyrLpState.data, gyrLp);

        s->restLastGyrLp.x = gyrLp[0];
        s->restLastGyrLp.y = gyrLp[1];
        s->restLastGyrLp.z = gyrLp[2];

        s->restLastSquaredDeviations.data[0] =
            (gyrVec[0] - gyrLp[0]) * (gyrVec[0] - gyrLp[0]) +
            (gyrVec[1] - gyrLp[1]) * (gyrVec[1] - gyrLp[1]) +
            (gyrVec[2] - gyrLp[2]) * (gyrVec[2] - gyrLp[2]);

        const float biasClip = p->biasClip * (constPI / 180.0f);
        const float restThGyr = p->restThGyr * (constPI / 180.0f);

        if (s->restLastSquaredDeviations.data[0] >= restThGyr * restThGyr ||
            fabsf(gyrLp[0]) > biasClip ||
            fabsf(gyrLp[1]) > biasClip ||
            fabsf(gyrLp[2]) > biasClip) {
            s->restT = 0.0f;
            s->restDetected = 0u;
        }
    }

    /* Remove estimated gyro bias */
    const axis3f_t gyrNoBias = {gyr.x - s->bias.x, gyr.y - s->bias.y, gyr.z - s->bias.z};
    const float gyrNoBiasVec[3] = {gyrNoBias.x, gyrNoBias.y, gyrNoBias.z};

    /* Gyroscope prediction step */
    const float gyrNorm = vqf_norm_f32(gyrNoBiasVec, 3u);
    const float angle = gyrNorm * c->gyrTs;

    if (gyrNorm > VQF_EPS) {
        const float half = 0.5f * angle;
        const float c0 = COS(half);
        const float s0 = SIN(half) / gyrNorm;

        quaternion_t step;
        step.q0 = c0;
        step.q1 = s0 * gyrNoBias.x;
        step.q2 = s0 * gyrNoBias.y;
        step.q3 = s0 * gyrNoBias.z;

        quaternionMult(&s->gyrQuat, &step, &s->gyrQuat);
        quaternionNorm(&s->gyrQuat);
    }
}

void AHRS_VQF_UpdateAcc(AHRS_VQF_t* vqf, axis3f_t acc)
{
    if (vqf == NULL) {
        return;
    }

    AHRS_VQF_Params_t* p = &vqf->params;
    AHRS_VQF_Coeffs_t* c = &vqf->coeffs;
    AHRS_VQF_State_t*  s = &vqf->state;

    if (acc.x == 0.0f && acc.y == 0.0f && acc.z == 0.0f) {
        return;
    }

    /* Convert body-NED sample to internal body-ENU */
    const axis3f_t accBody = vqf_body_ned_to_enu(acc);

    /* Rest detection (acc) */
    if (p->restBiasEstEnabled) {
        const float accVec[3] = {accBody.x, accBody.y, accBody.z};
        float accLp[3];

        vqf_filter_vec(accVec, 3u, p->restFilterTau, c->accTs,
                       c->restAccLpB.data, c->restAccLpA.data,
                       s->restAccLpState.data, accLp);

        s->restLastAccLp.x = accLp[0];
        s->restLastAccLp.y = accLp[1];
        s->restLastAccLp.z = accLp[2];

        s->restLastSquaredDeviations.data[1] =
            (accVec[0] - accLp[0]) * (accVec[0] - accLp[0]) +
            (accVec[1] - accLp[1]) * (accVec[1] - accLp[1]) +
            (accVec[2] - accLp[2]) * (accVec[2] - accLp[2]);

        if (s->restLastSquaredDeviations.data[1] >= p->restThAcc * p->restThAcc) {
            s->restT = 0.0f;
            s->restDetected = 0u;
        } else {
            s->restT += c->accTs;
            if (s->restT >= p->restMinT) {
                s->restDetected = 1u;
            }
        }
    }

    axis3f_t accEarth;

    /* Filter acc in inertial frame (earth-ENU) */
    vqf_quat_rotate(&s->gyrQuat, &accBody, &accEarth);

    float accEarthVec[3] = {accEarth.x, accEarth.y, accEarth.z};
    float accLpVec[3];

    vqf_filter_vec(accEarthVec, 3u, p->tauAcc, c->accTs,
                   c->accLpB.data, c->accLpA.data,
                   s->accLpState.data, accLpVec);

    s->lastAccLp.x = accLpVec[0];
    s->lastAccLp.y = accLpVec[1];
    s->lastAccLp.z = accLpVec[2];

    /* Transform to 6D earth frame and normalize */
    axis3f_t accEarth6D;
    vqf_quat_rotate(&s->accQuat, &s->lastAccLp, &accEarth6D);

    float accEarth6DVec[3] = {accEarth6D.x, accEarth6D.y, accEarth6D.z};
    vqf_normalize_f32(accEarth6DVec, 3u);

    accEarth6D.x = accEarth6DVec[0];
    accEarth6D.y = accEarth6DVec[1];
    accEarth6D.z = accEarth6DVec[2];

    /* Inclination correction */
    quaternion_t accCorrQuat;
    const float q_w = SQRT((accEarth6D.z + 1.0f) * 0.5f);

    if (q_w > 1e-6f) {
        accCorrQuat.q0 = q_w;
        accCorrQuat.q1 = 0.5f * accEarth6D.y / q_w;
        accCorrQuat.q2 = -0.5f * accEarth6D.x / q_w;
        accCorrQuat.q3 = 0.0f;
    } else {
        /* Near 180 degree correction */
        accCorrQuat.q0 = 0.0f;
        accCorrQuat.q1 = 1.0f;
        accCorrQuat.q2 = 0.0f;
        accCorrQuat.q3 = 0.0f;
    }

    quaternionMult(&accCorrQuat, &s->accQuat, &s->accQuat);
    quaternionNorm(&s->accQuat);

    /* Debug angular rate */
    s->lastAccCorrAngularRate = acosf(CONSTRAIN(accEarth6D.z, -1.0f, 1.0f)) / c->accTs;

    /* Bias estimation */
    if (p->motionBiasEstEnabled || p->restBiasEstEnabled) {
        const float biasClip = p->biasClip * (constPI / 180.0f);

        quaternion_t accGyrQuat;
        vqf_get_quat6d_enu(vqf, &accGyrQuat);

        float Rdata[9];
        matrix_t R = {Rdata, 3u, 3u};

        /* Rotation matrix R corresponding to accGyrQuat */
        const float qw = accGyrQuat.q0;
        const float qx = accGyrQuat.q1;
        const float qy = accGyrQuat.q2;
        const float qz = accGyrQuat.q3;

        Rdata[0] = 1.0f - 2.0f * (qy * qy) - 2.0f * (qz * qz);
        Rdata[1] = 2.0f * (qy * qx - qw * qz);
        Rdata[2] = 2.0f * (qw * qy + qz * qx);

        Rdata[3] = 2.0f * (qw * qz + qy * qx);
        Rdata[4] = 1.0f - 2.0f * (qx * qx) - 2.0f * (qz * qz);
        Rdata[5] = 2.0f * (qy * qz - qx * qw);

        Rdata[6] = 2.0f * (qz * qx - qw * qy);
        Rdata[7] = 2.0f * (qw * qx + qz * qy);
        Rdata[8] = 1.0f - 2.0f * (qx * qx) - 2.0f * (qy * qy);

        float biasLp[2];
        biasLp[0] = Rdata[0] * s->bias.x + Rdata[1] * s->bias.y + Rdata[2] * s->bias.z;
        biasLp[1] = Rdata[3] * s->bias.x + Rdata[4] * s->bias.y + Rdata[5] * s->bias.z;

        /* Low-pass filter R and R*b_hat */
        vqf_filter_vec(Rdata, 9u, p->tauAcc, c->accTs,
                       c->accLpB.data, c->accLpA.data,
                       s->motionBiasEstRLpState.data, Rdata);

        vqf_filter_vec(biasLp, 2u, p->tauAcc, c->accTs,
                       c->accLpB.data, c->accLpA.data,
                       s->motionBiasEstBiasLpState.data, biasLp);

        float w[3];
        float e[3];

        if (s->restDetected && p->restBiasEstEnabled) {
            e[0] = s->restLastGyrLp.x - s->bias.x;
            e[1] = s->restLastGyrLp.y - s->bias.y;
            e[2] = s->restLastGyrLp.z - s->bias.z;

            vqf_mat3_set_scaled_identity(1.0f, &R);

            w[0] = c->biasRestW;
            w[1] = c->biasRestW;
            w[2] = c->biasRestW;
        } else if (p->motionBiasEstEnabled) {
            /* Use normalized accEarth6D vector for motion update */
            e[0] = -accEarth6D.y / c->accTs + biasLp[0] - (Rdata[0] * s->bias.x + Rdata[1] * s->bias.y + Rdata[2] * s->bias.z);
            e[1] =  accEarth6D.x / c->accTs + biasLp[1] - (Rdata[3] * s->bias.x + Rdata[4] * s->bias.y + Rdata[5] * s->bias.z);
            e[2] = -(Rdata[6] * s->bias.x + Rdata[7] * s->bias.y + Rdata[8] * s->bias.z);

            w[0] = c->biasMotionW;
            w[1] = c->biasMotionW;
            w[2] = c->biasVerticalW;
        } else {
            w[0] = -1.0f; w[1] = -1.0f; w[2] = -1.0f;
        }

        /* Kalman update: step 1, P = P + V */
        if (s->biasP.data[0] < c->biasP0) {
            s->biasP.data[0] += c->biasV;
        }
        if (s->biasP.data[4] < c->biasP0) {
            s->biasP.data[4] += c->biasV;
        }
        if (s->biasP.data[8] < c->biasP0) {
            s->biasP.data[8] += c->biasV;
        }

        if (w[0] >= 0.0f) {
            /* Clip disagreement */
            e[0] = CONSTRAIN(e[0], -biasClip, biasClip);
            e[1] = CONSTRAIN(e[1], -biasClip, biasClip);
            e[2] = CONSTRAIN(e[2], -biasClip, biasClip);

            /* step 2: K = P R^T inv(W + R P R^T) */
            float K1Data[9];
            float SData[9];
            float invSData[9];
            float KData[9];
            float tmpData[9];
            float tmp2Data[9];

            matrix_t K1   = {K1Data,   3u, 3u};
            matrix_t S    = {SData,    3u, 3u};
            matrix_t invS = {invSData, 3u, 3u};
            matrix_t K    = {KData,    3u, 3u};
            matrix_t tmp  = {tmpData,  3u, 3u};
            matrix_t tmp2 = {tmp2Data, 3u, 3u};

            matrixMult_rhsT(&s->biasP, &R, &K1); /* K1 = P R^T */
            matrixMult(&R, &K1, &S);             /* S = R P R^T */
            S.data[0] += w[0];
            S.data[4] += w[1];
            S.data[8] += w[2];

            vqf_mat3_inv(&S, &invS);
            matrixMult(&K1, &invS, &K);          /* K = P R^T invS */

            /* step 3: bias = bias + K e */
            s->bias.x += K.data[0] * e[0] + K.data[1] * e[1] + K.data[2] * e[2];
            s->bias.y += K.data[3] * e[0] + K.data[4] * e[1] + K.data[5] * e[2];
            s->bias.z += K.data[6] * e[0] + K.data[7] * e[1] + K.data[8] * e[2];

            /* step 4: P = P - K R P */
            matrixMult(&K, &R, &tmp);            /* tmp = K R */
            matrixMult(&tmp, &s->biasP, &tmp2);  /* tmp2 = K R P */

            for (uint16_t i = 0; i < 9u; i++) {
                s->biasP.data[i] -= tmp2.data[i];
            }

            /* Clip bias estimate */
            s->bias.x = CONSTRAIN(s->bias.x, -biasClip, biasClip);
            s->bias.y = CONSTRAIN(s->bias.y, -biasClip, biasClip);
            s->bias.z = CONSTRAIN(s->bias.z, -biasClip, biasClip);
        }
    }
}

void AHRS_VQF_UpdateMag(AHRS_VQF_t* vqf, axis3f_t mag)
{
    if (vqf == NULL) {
        return;
    }

    AHRS_VQF_Params_t* p = &vqf->params;
    AHRS_VQF_Coeffs_t* c = &vqf->coeffs;
    AHRS_VQF_State_t*  s = &vqf->state;

    if (mag.x == 0.0f && mag.y == 0.0f && mag.z == 0.0f) {
        return;
    }

    /* Convert body-NED sample to internal body-ENU */
    const axis3f_t magBody = vqf_body_ned_to_enu(mag);

    axis3f_t magEarth;

    /* Bring magnetometer measurement into 6D earth frame */
    quaternion_t accGyrQuat;
    vqf_get_quat6d_enu(vqf, &accGyrQuat);
    vqf_quat_rotate(&accGyrQuat, &magBody, &magEarth);

    if (p->magDistRejectionEnabled) {
        const float magEarthVec[3] = {magEarth.x, magEarth.y, magEarth.z};
        s->magNormDip.data[0] = vqf_norm_f32(magEarthVec, 3u);

        if (s->magNormDip.data[0] > VQF_EPS) {
            s->magNormDip.data[1] = -asinf(CONSTRAIN(magEarth.z / s->magNormDip.data[0], -1.0f, 1.0f));
        } else {
            s->magNormDip.data[1] = 0.0f;
        }

        if (p->magCurrentTau > 0.0f) {
            vqf_filter_vec(s->magNormDip.data, 2u, p->magCurrentTau, c->magTs,
                           c->magNormDipLpB.data, c->magNormDipLpA.data,
                           s->magNormDipLpState.data, s->magNormDip.data);
        }

        /* Magnetic disturbance detection */
        if (fabsf(s->magNormDip.data[0] - s->magRefNorm) < p->magNormTh * s->magRefNorm &&
            fabsf(s->magNormDip.data[1] - s->magRefDip) < p->magDipTh * (constPI / 180.0f)) {
            s->magUndisturbedT += c->magTs;
            if (s->magUndisturbedT >= p->magMinUndisturbedTime) {
                s->magDistDetected = 0u;
                s->magRefNorm += c->kMagRef * (s->magNormDip.data[0] - s->magRefNorm);
                s->magRefDip  += c->kMagRef * (s->magNormDip.data[1] - s->magRefDip);
            }
        } else {
            s->magUndisturbedT = 0.0f;
            s->magDistDetected = 1u;
        }

        /* New magnetic field acceptance */
        if (fabsf(s->magNormDip.data[0] - s->magCandidateNorm) < p->magNormTh * s->magCandidateNorm &&
            fabsf(s->magNormDip.data[1] - s->magCandidateDip) < p->magDipTh * (constPI / 180.0f)) {

            const float restLastGyrLpVec[3] = {s->restLastGyrLp.x, s->restLastGyrLp.y, s->restLastGyrLp.z};
            if (vqf_norm_f32(restLastGyrLpVec, 3u) >= p->magNewMinGyr * (constPI / 180.0f)) {
                s->magCandidateT += c->magTs;
            }

            s->magCandidateNorm += c->kMagRef * (s->magNormDip.data[0] - s->magCandidateNorm);
            s->magCandidateDip  += c->kMagRef * (s->magNormDip.data[1] - s->magCandidateDip);

            if (s->magDistDetected &&
                (s->magCandidateT >= p->magNewTime ||
                 (s->magRefNorm == 0.0f && s->magCandidateT >= p->magNewFirstTime))) {

                s->magRefNorm = s->magCandidateNorm;
                s->magRefDip = s->magCandidateDip;
                s->magDistDetected = 0u;
                s->magUndisturbedT = p->magMinUndisturbedTime;
            }
        } else {
            s->magCandidateT = 0.0f;
            s->magCandidateNorm = s->magNormDip.data[0];
            s->magCandidateDip = s->magNormDip.data[1];
        }
    }

    /* Disagreement angle */
    s->lastMagDisAngle = atan2f(magEarth.x, magEarth.y) - s->delta;

    if (s->lastMagDisAngle > constPI) {
        s->lastMagDisAngle -= 2.0f * constPI;
    } else if (s->lastMagDisAngle < -constPI) {
        s->lastMagDisAngle += 2.0f * constPI;
    }

    float k = c->kMag;

    if (p->magDistRejectionEnabled) {
        if (s->magDistDetected) {
            if (s->magRejectT <= p->magMaxRejectionTime) {
                s->magRejectT += c->magTs;
                k = 0.0f;
            } else {
                k /= p->magRejectionFactor;
            }
        } else {
            s->magRejectT = MAX(s->magRejectT - p->magRejectionFactor * c->magTs, 0.0f);
        }
    }

    /* Fast initial convergence */
    if (s->kMagInit != 0.0f) {
        if (k < s->kMagInit) {
            k = s->kMagInit;
        }

        s->kMagInit = s->kMagInit / (s->kMagInit + 1.0f);

        if (s->kMagInit * p->tauMag < c->magTs) {
            s->kMagInit = 0.0f;
        }
    }

    /* First-order filter step */
    s->delta += k * s->lastMagDisAngle;
    s->lastMagCorrAngularRate = k * s->lastMagDisAngle / c->magTs;

    if (s->delta > constPI) {
        s->delta -= 2.0f * constPI;
    } else if (s->delta < -constPI) {
        s->delta += 2.0f * constPI;
    }
}

/* Outputs -------------------------------------------------------------------*/

void AHRS_VQF_GetQuat3D(const AHRS_VQF_t* vqf, quaternion_t* q_out)
{
    if (vqf == NULL || q_out == NULL) {
        return;
    }

    quaternion_t q_enu;
    vqf_get_quat3d_enu(vqf, &q_enu);
    vqf_quat_enu_to_ned(&q_enu, q_out);
}

void AHRS_VQF_GetQuat6D(const AHRS_VQF_t* vqf, quaternion_t* q_out)
{
    if (vqf == NULL || q_out == NULL) {
        return;
    }

    quaternion_t q_enu;
    vqf_get_quat6d_enu(vqf, &q_enu);
    vqf_quat_enu_to_ned(&q_enu, q_out);
}

void AHRS_VQF_GetQuat9D(const AHRS_VQF_t* vqf, quaternion_t* q_out)
{
    if (vqf == NULL || q_out == NULL) {
        return;
    }

    quaternion_t q_enu;
    vqf_get_quat9d_enu(vqf, &q_enu);
    vqf_quat_enu_to_ned(&q_enu, q_out);
}

float AHRS_VQF_GetDelta(const AHRS_VQF_t* vqf)
{
    if (vqf == NULL) {
        return 0.0f;
    }

    /* Convert ENU +Z (Up) to NED +Z (Down) */
    return -vqf->state.delta;
}

float AHRS_VQF_GetBiasEstimate(const AHRS_VQF_t* vqf, axis3f_t* bias_out)
{
    if (vqf == NULL) {
        return 0.0f;
    }

    const AHRS_VQF_State_t* s = &vqf->state;

    if (bias_out != NULL) {
        /* Convert internal body-ENU bias to body-NED */
        *bias_out = vqf_body_enu_to_ned(s->bias);
    }

    const float p0 = s->biasP.data[0];
    const float p1 = s->biasP.data[4];
    const float p2 = s->biasP.data[8];
    const float pMax = MAX(p0, MAX(p1, p2));

    /* Convert from internal 0.01deg/s scaling to rad/s (as in the reference implementation) */
    return SQRT(pMax) * (constPI / (100.0f * 180.0f));
}

uint8_t AHRS_VQF_GetRestDetected(const AHRS_VQF_t* vqf)
{
    if (vqf == NULL) {
        return 0u;
    }
    return vqf->state.restDetected;
}

uint8_t AHRS_VQF_GetMagDistDetected(const AHRS_VQF_t* vqf)
{
    if (vqf == NULL) {
        return 0u;
    }
    return vqf->state.magDistDetected;
}

void AHRS_VQF_GetRelativeRestDeviations(const AHRS_VQF_t* vqf, matrix_t* out)
{
    if (vqf == NULL || out == NULL || out->rows != 2u || out->cols != 1u || out->data == NULL) {
        return;
    }

    const float thGyr = vqf->params.restThGyr * (constPI / 180.0f);
    const float thAcc = vqf->params.restThAcc;

    out->data[0] = (thGyr > 0.0f) ? (SQRT(vqf->state.restLastSquaredDeviations.data[0]) / thGyr) : 0.0f;
    out->data[1] = (thAcc > 0.0f) ? (SQRT(vqf->state.restLastSquaredDeviations.data[1]) / thAcc) : 0.0f;
}

float AHRS_VQF_GetMagRefNorm(const AHRS_VQF_t* vqf)
{
    if (vqf == NULL) {
        return 0.0f;
    }
    return vqf->state.magRefNorm;
}

float AHRS_VQF_GetMagRefDip(const AHRS_VQF_t* vqf)
{
    if (vqf == NULL) {
        return 0.0f;
    }
    return vqf->state.magRefDip;
}

/* Setters -------------------------------------------------------------------*/

void AHRS_VQF_SetTauAcc(AHRS_VQF_t* vqf, float tauAcc_s)
{
    if (vqf == NULL) {
        return;
    }

    AHRS_VQF_Params_t* p = &vqf->params;

    if (p->tauAcc == tauAcc_s) {
        return;
    }

    float b_new[3];
    float a_new[2];

    vqf_filter_coeffs(tauAcc_s, vqf->coeffs.accTs, b_new, a_new);

    {
        float lastAccLp[3] = {vqf->state.lastAccLp.x, vqf->state.lastAccLp.y, vqf->state.lastAccLp.z};

        vqf_filter_adapt_state_for_coeff_change(lastAccLp, 3u,
                                                vqf->coeffs.accLpB.data, vqf->coeffs.accLpA.data,
                                                b_new, a_new,
                                                vqf->state.accLpState.data);
    }

    if (vqf->params.motionBiasEstEnabled) {
        float r_last[9];
        float bias_last[2];

        for (uint16_t i = 0; i < 9u; i++) {
            r_last[i] = vqf->state.motionBiasEstRLpState.data[2u * i];
        }
        for (uint16_t i = 0; i < 2u; i++) {
            bias_last[i] = vqf->state.motionBiasEstBiasLpState.data[2u * i];
        }

        vqf_filter_adapt_state_for_coeff_change(r_last, 9u,
                                                vqf->coeffs.accLpB.data, vqf->coeffs.accLpA.data,
                                                b_new, a_new,
                                                vqf->state.motionBiasEstRLpState.data);

        vqf_filter_adapt_state_for_coeff_change(bias_last, 2u,
                                                vqf->coeffs.accLpB.data, vqf->coeffs.accLpA.data,
                                                b_new, a_new,
                                                vqf->state.motionBiasEstBiasLpState.data);
    }

    memcpy(vqf->coeffs.accLpB.data, b_new, sizeof(b_new));
    memcpy(vqf->coeffs.accLpA.data, a_new, sizeof(a_new));

    p->tauAcc = tauAcc_s;
}

void AHRS_VQF_SetTauMag(AHRS_VQF_t* vqf, float tauMag_s)
{
    if (vqf == NULL) {
        return;
    }

    AHRS_VQF_Params_t* p = &vqf->params;

    if (p->tauMag == tauMag_s) {
        return;
    }

    vqf->coeffs.kMag = vqf_gain_from_tau(tauMag_s, vqf->coeffs.magTs);
    p->tauMag = tauMag_s;
}

void AHRS_VQF_SetMotionBiasEstEnabled(AHRS_VQF_t* vqf, uint8_t enabled)
{
    if (vqf == NULL) {
        return;
    }

    vqf->params.motionBiasEstEnabled = enabled ? 1u : 0u;

    for (uint16_t i = 0; i < 18u; i++) {
        vqf->state.motionBiasEstRLpState.data[i] = VQF_NAN;
    }
    for (uint16_t i = 0; i < 4u; i++) {
        vqf->state.motionBiasEstBiasLpState.data[i] = VQF_NAN;
    }
}

void AHRS_VQF_SetRestBiasEstEnabled(AHRS_VQF_t* vqf, uint8_t enabled)
{
    if (vqf == NULL) {
        return;
    }

    vqf->params.restBiasEstEnabled = enabled ? 1u : 0u;

    vqf->state.restDetected = 0u;
    vqf->state.restT = 0.0f;

    vqf->state.restLastSquaredDeviations.data[0] = 0.0f;
    vqf->state.restLastSquaredDeviations.data[1] = 0.0f;

    for (uint16_t i = 0; i < 6u; i++) {
        vqf->state.restGyrLpState.data[i] = VQF_NAN;
        vqf->state.restAccLpState.data[i] = VQF_NAN;
    }
}

void AHRS_VQF_SetMagDistRejectionEnabled(AHRS_VQF_t* vqf, uint8_t enabled)
{
    if (vqf == NULL) {
        return;
    }

    vqf->params.magDistRejectionEnabled = enabled ? 1u : 0u;

    vqf->state.magDistDetected = 1u;
    vqf->state.magUndisturbedT = 0.0f;
    vqf->state.magRejectT = vqf->params.magMaxRejectionTime;

    vqf->state.magCandidateNorm = -1.0f;
    vqf->state.magCandidateDip = 0.0f;
    vqf->state.magCandidateT = 0.0f;

    for (uint16_t i = 0; i < 4u; i++) {
        vqf->state.magNormDipLpState.data[i] = VQF_NAN;
    }
}

void AHRS_VQF_SetRestDetectionThresholds(AHRS_VQF_t* vqf, float thGyr, float thAcc)
{
    if (vqf == NULL) {
        return;
    }

    vqf->params.restThGyr = thGyr;
    vqf->params.restThAcc = thAcc;
}

void AHRS_VQF_SetMagRef(AHRS_VQF_t* vqf, float norm, float dip)
{
    if (vqf == NULL) {
        return;
    }

    vqf->state.magRefNorm = norm;
    vqf->state.magRefDip = dip;
}

void AHRS_VQF_SetBiasEstimate(AHRS_VQF_t* vqf, axis3f_t bias, float sigma)
{
    if (vqf == NULL) {
        return;
    }

    /* Convert provided body-NED bias to internal body-ENU */
    vqf->state.bias = vqf_body_ned_to_enu(bias);

    if (sigma > 0.0f) {
        const float tmp = sigma * (180.0f * 100.0f / constPI); /* rad/s -> 0.01deg/s */
        const float p0 = tmp * tmp;
        vqf_mat3_set_scaled_identity(p0, &vqf->state.biasP);
    }
}
