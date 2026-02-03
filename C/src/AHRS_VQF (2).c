// Copyright (c) 2024 Hugo Chiang
// SPDX-License-Identifier: MIT
//
// NOTE: This file is a refactor of the original vqf.c to match the InertialEstimators
// C library structure (instance-based state, consistent naming).
//
// Coordinate convention
// - Public API expects sensor samples in body-NED (X North/Forward, Y East/Right, Z Down)
//   when configAHRS_VQF_INPUT_BODY_NED=1 (default).
// - Internally VQF runs in an ENU (+Z Up) convention; all inputs/outputs are mapped so
//   the external behaviour is body->NED.

#include "AHRS_VQF.h"

#include "basicMath.h"
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define VQF_PI_F   (3.14159265358979323846f)
#define VQF_SQRT2  (1.41421356237309504880f)
#define VQF_EPS    (FLT_EPSILON)
#define VQF_NAN    (NAN)

#if defined(USE_FAST_MATH)
#include "basicMath.h"
static inline float vqf_sin(float x) { return fastSin(x); }
static inline float vqf_cos(float x) { return fastCos(x); }
static inline float vqf_sqrt(float x) { return fastSqrt(x); }
static inline float vqf_inv_sqrt(float x) { return fastInvSqrt(x); }
#else
static inline float vqf_sin(float x) { return sinf(x); }
static inline float vqf_cos(float x) { return cosf(x); }
static inline float vqf_sqrt(float x) { return sqrtf(x); }
static inline float vqf_inv_sqrt(float x) { return 1.0f / sqrtf(x); }
#endif

static inline float vqf_tan(float x) { return tanf(x); }
static inline float vqf_exp(float x) { return expf(x); }
static inline float vqf_atan2(float y, float x) { return atan2f(y, x); }
static inline float vqf_asin(float x) { return asinf(x); }
static inline float vqf_acos(float x) { return acosf(x); }
static inline float vqf_fabs(float x) { return fabsf(x); }

// -----------------------------------------------------------------------------
// Small utilities
// -----------------------------------------------------------------------------

static void vqf_fill_f32(float* dst, size_t n, float val)
{
    for (size_t i = 0; i < n; i++) {
        dst[i] = val;
    }
}

static void vqf_fill_f64(double* dst, size_t n, double val)
{
    for (size_t i = 0; i < n; i++) {
        dst[i] = val;
    }
}

static inline float vqf_square(float x) { return x * x; }

static inline float vqf_min(float a, float b) { return (a < b) ? a : b; }
static inline float vqf_max(float a, float b) { return (a > b) ? a : b; }

static float vqf_norm(const float* vec, size_t n)
{
    float s = 0.0f;
    for (size_t i = 0; i < n; i++) {
        s += vec[i] * vec[i];
    }
    return vqf_sqrt(s);
}

static void vqf_normalize(float* vec, size_t n)
{
    float nrm = vqf_norm(vec, n);
    if (nrm <= VQF_EPS) {
        return;
    }
    float inv = 1.0f / nrm;
    for (size_t i = 0; i < n; i++) {
        vec[i] *= inv;
    }
}

static void vqf_clip(float* vec, size_t n, float minVal, float maxVal)
{
    for (size_t i = 0; i < n; i++) {
        if (vec[i] < minVal) {
            vec[i] = minVal;
        } else if (vec[i] > maxVal) {
            vec[i] = maxVal;
        }
    }
}

// -----------------------------------------------------------------------------
// Quaternion math (scalar-first [w x y z])
// -----------------------------------------------------------------------------

static void vqf_quat_set_identity(float q[4])
{
    q[0] = 1.0f;
    q[1] = 0.0f;
    q[2] = 0.0f;
    q[3] = 0.0f;
}


static void vqf_quat_multiply(const float q1[4], const float q2[4], float out[4])
{
    const float w = q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2] - q1[3] * q2[3];
    const float x = q1[0] * q2[1] + q1[1] * q2[0] + q1[2] * q2[3] - q1[3] * q2[2];
    const float y = q1[0] * q2[2] - q1[1] * q2[3] + q1[2] * q2[0] + q1[3] * q2[1];
    const float z = q1[0] * q2[3] + q1[1] * q2[2] - q1[2] * q2[1] + q1[3] * q2[0];

    out[0] = w;
    out[1] = x;
    out[2] = y;
    out[3] = z;
}

static void vqf_quat_rotate(const float q[4], const float v[3], float out[3])
{
    const float x = (1.0f - 2.0f * q[2] * q[2] - 2.0f * q[3] * q[3]) * v[0]
                  + 2.0f * v[1] * (q[2] * q[1] - q[0] * q[3])
                  + 2.0f * v[2] * (q[0] * q[2] + q[3] * q[1]);
    const float y = 2.0f * v[0] * (q[0] * q[3] + q[2] * q[1])
                  + (1.0f - 2.0f * q[1] * q[1] - 2.0f * q[3] * q[3]) * v[1]
                  + 2.0f * v[2] * (q[2] * q[3] - q[1] * q[0]);
    const float z = 2.0f * v[0] * (q[3] * q[1] - q[0] * q[2])
                  + 2.0f * v[1] * (q[0] * q[1] + q[3] * q[2])
                  + (1.0f - 2.0f * q[1] * q[1] - 2.0f * q[2] * q[2]) * v[2];

    out[0] = x;
    out[1] = y;
    out[2] = z;
}

static void vqf_quat_apply_delta(const float q[4], float delta, float out[4])
{
    // out = [cos(delta/2), 0, 0, sin(delta/2)] ⊗ q
    const float c = vqf_cos(delta * 0.5f);
    const float s = vqf_sin(delta * 0.5f);

    const float w = c * q[0] - s * q[3];
    const float x = c * q[1] - s * q[2];
    const float y = c * q[2] + s * q[1];
    const float z = c * q[3] + s * q[0];

    out[0] = w;
    out[1] = x;
    out[2] = y;
    out[3] = z;
}

// ENU <-> NED axis mapping (swap X/Y and flip Z):
// v_ned = P * v_enu,  P = [[0,1,0],[1,0,0],[0,0,-1]]
//
// If inputs are body-NED (configAHRS_VQF_INPUT_BODY_NED=1), internal VQF runs in ENU
// and the output quaternion must convert both frames:
//   q_bodyNED_to_navNED = qP ⊗ q_bodyENU_to_navENU ⊗ qP
//
// Otherwise (body already ENU), only the navigation frame conversion is applied:
//   q_bodyENU_to_navNED = qP ⊗ q_bodyENU_to_navENU

static inline void vqf_vec_body_to_enu(const float v_body[3], float v_enu[3])
{
#if configAHRS_VQF_INPUT_BODY_NED
    // body-NED -> body-ENU: [E, N, U] = [Y, X, -Z]
    v_enu[0] = v_body[1];
    v_enu[1] = v_body[0];
    v_enu[2] = -v_body[2];
#else
    // body already expressed as ENU
    v_enu[0] = v_body[0];
    v_enu[1] = v_body[1];
    v_enu[2] = v_body[2];
#endif
}

static inline void vqf_vec_enu_to_body(const float v_enu[3], float v_body[3])
{
#if configAHRS_VQF_INPUT_BODY_NED
    // body-ENU -> body-NED: [N, E, D] = [Y, X, -Z]
    v_body[0] = v_enu[1];
    v_body[1] = v_enu[0];
    v_body[2] = -v_enu[2];
#else
    v_body[0] = v_enu[0];
    v_body[1] = v_enu[1];
    v_body[2] = v_enu[2];
#endif
}

static void vqf_quat_enu_to_ned(const float q_enu[4], float q_ned[4])
{
    const float s = 0.70710678118654752440f;
    const float qP[4] = {0.0f, s, s, 0.0f};

#if configAHRS_VQF_INPUT_BODY_NED
    // q_ned = qP ⊗ q_enu ⊗ qP  (convert both body and navigation frames)
    float tmp[4];
    vqf_quat_multiply(qP, q_enu, tmp);
    vqf_quat_multiply(tmp, qP, q_ned);
#else
    // q_ned = qP ⊗ q_enu  (convert navigation frame only)
    vqf_quat_multiply(qP, q_enu, q_ned);
#endif

    // numerical noise: keep unit length
    vqf_normalize(q_ned, 4);
}

// -----------------------------------------------------------------------------
// 2nd order Butterworth LP filter helpers (per original VQF)
// -----------------------------------------------------------------------------

static float vqf_gain_from_tau(float tau, float Ts)
{
    if (tau < 0.0f) {
        return 0.0f; // k=0 for negative tau (disable update)
    } else if (tau == 0.0f) {
        return 1.0f; // k=1 for tau=0
    } else {
        return 1.0f - vqf_exp(-Ts / tau); // fc = 1/(2*pi*tau)
    }
}

static void vqf_filter_coeffs(float tau, float Ts, double outB[3], double outA[2])
{
    // second order Butterworth filter based on https://stackoverflow.com/a/52764064
    // time constant of damped, non-oscillating part of step response
    const double fc = (VQF_SQRT2 / (2.0 * VQF_PI_F)) / (double)tau;

    const double C = (double)vqf_tan((float)(VQF_PI_F * (float)fc * Ts));
    const double D = C * C + sqrt(2.0) * C + 1.0;

    const double b0 = (C * C) / D;
    outB[0] = b0;
    outB[1] = 2.0 * b0;
    outB[2] = b0;

    // a0 = 1.0
    outA[0] = 2.0 * (C * C - 1.0) / D;                 // a1
    outA[1] = (C * C - sqrt(2.0) * C + 1.0) / D;       // a2
}

static void vqf_filter_initial_state(float x0, const double b[3], const double a[2], double out[2])
{
    // Initial state for steady state (equivalent to scipy.signal.lfilter_zi).
    out[0] = (double)x0 * (1.0 - b[0]);
    out[1] = (double)x0 * (b[2] - a[1]);
}

static void vqf_filter_adapt_state_for_coeff_change(
    const float last_y[], size_t n,
    const double b_old[3], const double a_old[2],
    const double b_new[3], const double a_new[2],
    double state[])
{
    if (isnan(state[0])) {
        return;
    }
    for (size_t i = 0; i < n; i++) {
        state[0 + 2 * i] = state[0 + 2 * i] + (b_old[0] - b_new[0]) * (double)last_y[i];
        state[1 + 2 * i] = state[1 + 2 * i] + (b_old[1] - b_new[1] - a_old[0] + a_new[0]) * (double)last_y[i];
    }
}

static float vqf_filter_step(float x, const double b[3], const double a[2], double state[2])
{
    // difference equations based on scipy.signal.lfilter documentation
    // assumes that a0 == 1.0
    const double y = b[0] * (double)x + state[0];
    state[0] = b[1] * (double)x - a[0] * y + state[1];
    state[1] = b[2] * (double)x - a[1] * y;
    return (float)y;
}

static void vqf_filter_vec(
    const float x[], size_t n, float tau, float Ts,
    const double b[3], const double a[2], double state[], float out[])
{
    // assert(n >= 2);

    // Initialization phase:
    // To avoid depending on a single sample, average the first samples (for duration tau)
    // and then use this average to calculate the filter initial state.
    if (isnan(state[0])) {
        if (isnan(state[1])) { // first sample
            state[1] = 0.0; // state[1] stores sample count
            for (size_t i = 0; i < n; i++) {
                state[2 + i] = 0.0; // state[2+i] stores sum
            }
        }

        state[1] += 1.0;
        for (size_t i = 0; i < n; i++) {
            state[2 + i] += (double)x[i];
            out[i] = (float)(state[2 + i] / state[1]);
        }

        if ((float)(state[1]) * Ts >= tau) {
            for (size_t i = 0; i < n; i++) {
                vqf_filter_initial_state(out[i], b, a, state + 2 * i);
            }
        }
        return;
    }

    for (size_t i = 0; i < n; i++) {
        out[i] = vqf_filter_step(x[i], b, a, state + 2 * i);
    }
}
// -----------------------------------------------------------------------------
// 3x3 matrix helpers
// -----------------------------------------------------------------------------

static void vqf_mat3_set_scaled_identity(float s, float out[9])
{
    vqf_fill_f32(out, 9, 0.0f);
    out[0] = s;
    out[4] = s;
    out[8] = s;
}

static void vqf_mat3_multiply(const float A[9], const float B[9], float out[9])
{
    out[0] = A[0] * B[0] + A[1] * B[3] + A[2] * B[6];
    out[1] = A[0] * B[1] + A[1] * B[4] + A[2] * B[7];
    out[2] = A[0] * B[2] + A[1] * B[5] + A[2] * B[8];

    out[3] = A[3] * B[0] + A[4] * B[3] + A[5] * B[6];
    out[4] = A[3] * B[1] + A[4] * B[4] + A[5] * B[7];
    out[5] = A[3] * B[2] + A[4] * B[5] + A[5] * B[8];

    out[6] = A[6] * B[0] + A[7] * B[3] + A[8] * B[6];
    out[7] = A[6] * B[1] + A[7] * B[4] + A[8] * B[7];
    out[8] = A[6] * B[2] + A[7] * B[5] + A[8] * B[8];
}

static void vqf_mat3_multiply_tps_first(const float A[9], const float B[9], float out[9])
{
    // out = A^T * B
    out[0] = A[0] * B[0] + A[3] * B[3] + A[6] * B[6];
    out[1] = A[0] * B[1] + A[3] * B[4] + A[6] * B[7];
    out[2] = A[0] * B[2] + A[3] * B[5] + A[6] * B[8];

    out[3] = A[1] * B[0] + A[4] * B[3] + A[7] * B[6];
    out[4] = A[1] * B[1] + A[4] * B[4] + A[7] * B[7];
    out[5] = A[1] * B[2] + A[4] * B[5] + A[7] * B[8];

    out[6] = A[2] * B[0] + A[5] * B[3] + A[8] * B[6];
    out[7] = A[2] * B[1] + A[5] * B[4] + A[8] * B[7];
    out[8] = A[2] * B[2] + A[5] * B[5] + A[8] * B[8];
}

static void vqf_mat3_multiply_tps_second(const float A[9], const float B[9], float out[9])
{
    // out = A * B^T
    out[0] = A[0] * B[0] + A[1] * B[1] + A[2] * B[2];
    out[1] = A[0] * B[3] + A[1] * B[4] + A[2] * B[5];
    out[2] = A[0] * B[6] + A[1] * B[7] + A[2] * B[8];

    out[3] = A[3] * B[0] + A[4] * B[1] + A[5] * B[2];
    out[4] = A[3] * B[3] + A[4] * B[4] + A[5] * B[5];
    out[5] = A[3] * B[6] + A[4] * B[7] + A[5] * B[8];

    out[6] = A[6] * B[0] + A[7] * B[1] + A[8] * B[2];
    out[7] = A[6] * B[3] + A[7] * B[4] + A[8] * B[5];
    out[8] = A[6] * B[6] + A[7] * B[7] + A[8] * B[8];
}

static uint8_t vqf_mat3_inv(const float A[9], float out[9])
{
    const float det =
        A[0] * (A[4] * A[8] - A[5] * A[7]) -
        A[1] * (A[3] * A[8] - A[5] * A[6]) +
        A[2] * (A[3] * A[7] - A[4] * A[6]);

    if (vqf_fabs(det) < 1e-12f) {
        return 0;
    }

    const float invDet = 1.0f / det;

    out[0] =  (A[4] * A[8] - A[5] * A[7]) * invDet;
    out[1] = -(A[1] * A[8] - A[2] * A[7]) * invDet;
    out[2] =  (A[1] * A[5] - A[2] * A[4]) * invDet;

    out[3] = -(A[3] * A[8] - A[5] * A[6]) * invDet;
    out[4] =  (A[0] * A[8] - A[2] * A[6]) * invDet;
    out[5] = -(A[0] * A[5] - A[2] * A[3]) * invDet;

    out[6] =  (A[3] * A[7] - A[4] * A[6]) * invDet;
    out[7] = -(A[0] * A[7] - A[1] * A[6]) * invDet;
    out[8] =  (A[0] * A[4] - A[1] * A[3]) * invDet;

    return 1;
}

// -----------------------------------------------------------------------------
// Internal setup / defaults
// -----------------------------------------------------------------------------

static void vqf_init_params(AHRS_VQF_Params_t* p)
{
    p->tauAcc = 3.0f;
    p->tauMag = 9.0f;
    p->motionBiasEstEnabled = 1;
    p->restBiasEstEnabled = 1;
    p->magDistRejectionEnabled = 1;
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

    vqf_quat_set_identity(s->gyrQuat);
    vqf_quat_set_identity(s->accQuat);
    s->delta = 0.0f;

    s->restDetected = 0;
    s->magDistDetected = 1;

    vqf_fill_f32(s->lastAccLp, 3, 0.0f);
    vqf_fill_f64(s->accLpState, 3 * 2, (double)VQF_NAN);
    s->lastAccCorrAngularRate = 0.0f;

    s->kMagInit = 1.0f;
    s->lastMagDisAngle = 0.0f;
    s->lastMagCorrAngularRate = 0.0f;

    vqf_fill_f32(s->bias, 3, 0.0f);
    vqf_mat3_set_scaled_identity(vqf->coeffs.biasP0, s->biasP);

    vqf_fill_f64(s->motionBiasEstRLpState, 9 * 2, (double)VQF_NAN);
    vqf_fill_f64(s->motionBiasEstBiasLpState, 2 * 2, (double)VQF_NAN);

    vqf_fill_f32(s->restLastSquaredDeviations, 2, 0.0f);
    s->restT = 0.0f;
    vqf_fill_f32(s->restLastGyrLp, 3, 0.0f);
    vqf_fill_f64(s->restGyrLpState, 3 * 2, (double)VQF_NAN);
    vqf_fill_f32(s->restLastAccLp, 3, 0.0f);
    vqf_fill_f64(s->restAccLpState, 3 * 2, (double)VQF_NAN);

    s->magRefNorm = 0.0f;
    s->magRefDip = 0.0f;
    s->magUndisturbedT = 0.0f;
    s->magRejectT = vqf->params.magMaxRejectionTime;

    s->magCandidateNorm = -1.0f;
    s->magCandidateDip = 0.0f;
    s->magCandidateT = 0.0f;

    vqf_fill_f32(s->magNormDip, 2, 0.0f);
    vqf_fill_f64(s->magNormDipLpState, 2 * 2, (double)VQF_NAN);
}

static void vqf_setup(AHRS_VQF_t* vqf)
{
    AHRS_VQF_Params_t* p = &vqf->params;
    AHRS_VQF_Coeffs_t* c = &vqf->coeffs;

    // Core accelerometer LP filter
    vqf_filter_coeffs(p->tauAcc, c->accTs, c->accLpB, c->accLpA);

    // Magnetic update gain
    c->kMag = vqf_gain_from_tau(p->tauMag, c->magTs);

    // Gyro bias estimator
    c->biasP0 = vqf_square(p->biasSigmaInit * 100.0f);
    c->biasV  = vqf_square(0.1f * 100.0f) * c->accTs / p->biasForgettingTime;

    const float pMotion = vqf_square(p->biasSigmaMotion * 100.0f);
    c->biasMotionW = vqf_square(pMotion) / c->biasV + pMotion;
    c->biasVerticalW = c->biasMotionW / vqf_max(p->biasVerticalForgettingFactor, 1e-10f);

    const float pRest = vqf_square(p->biasSigmaRest * 100.0f);
    c->biasRestW = vqf_square(pRest) / c->biasV + pRest;

    // Rest detection filters
    vqf_filter_coeffs(p->restFilterTau, c->gyrTs, c->restGyrLpB, c->restGyrLpA);
    vqf_filter_coeffs(p->restFilterTau, c->accTs, c->restAccLpB, c->restAccLpA);

    // Magnetic reference tracking
    c->kMagRef = vqf_gain_from_tau(p->magRefTau, c->magTs);

    if (p->magCurrentTau > 0.0f) {
        vqf_filter_coeffs(p->magCurrentTau, c->magTs, c->magNormDipLpB, c->magNormDipLpA);
    } else {
        vqf_fill_f64(c->magNormDipLpB, 3, (double)VQF_NAN);
        vqf_fill_f64(c->magNormDipLpA, 2, (double)VQF_NAN);
    }

    vqf_reset_state(vqf);
}

// -----------------------------------------------------------------------------
// Public API
// -----------------------------------------------------------------------------

void AHRS_VQF_Init(AHRS_VQF_t* vqf, float gyrTs_s, float accTs_s, float magTs_s)
{
    if (vqf == NULL) {
        return;
    }

    memset(vqf, 0, sizeof(*vqf));
    vqf_init_params(&vqf->params);

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

// -----------------------------------------------------------------------------
// Internal getters (ENU). Public getters convert to NED.
// -----------------------------------------------------------------------------

static void vqf_get_quat3d_enu(const AHRS_VQF_t* vqf, float out[4])
{
    memcpy(out, vqf->state.gyrQuat, sizeof(vqf->state.gyrQuat));
}

static void vqf_get_quat6d_enu(const AHRS_VQF_t* vqf, float out[4])
{
    vqf_quat_multiply(vqf->state.accQuat, vqf->state.gyrQuat, out);
}

static void vqf_get_quat9d_enu(const AHRS_VQF_t* vqf, float out[4])
{
    float q6d[4];
    vqf_get_quat6d_enu(vqf, q6d);
    vqf_quat_apply_delta(q6d, vqf->state.delta, out);
}

void AHRS_VQF_GetQuat3D(const AHRS_VQF_t* vqf, float q_out[4])
{
    float q_enu[4];
    vqf_get_quat3d_enu(vqf, q_enu);
    vqf_quat_enu_to_ned(q_enu, q_out);
}

void AHRS_VQF_GetQuat6D(const AHRS_VQF_t* vqf, float q_out[4])
{
    float q_enu[4];
    vqf_get_quat6d_enu(vqf, q_enu);
    vqf_quat_enu_to_ned(q_enu, q_out);
}

void AHRS_VQF_GetQuat9D(const AHRS_VQF_t* vqf, float q_out[4])
{
    float q_enu[4];
    vqf_get_quat9d_enu(vqf, q_enu);
    vqf_quat_enu_to_ned(q_enu, q_out);
}

float AHRS_VQF_GetDelta(const AHRS_VQF_t* vqf)
{
    if (vqf == NULL) {
        return 0.0f;
    }

    // Internal delta is defined around ENU +Z (Up). In NED, +Z is Down.
    // Represent the same physical yaw correction as rotation about +Down:
    return -vqf->state.delta;
}

// -----------------------------------------------------------------------------
// Parameter setters (kept close to original behavior)
// -----------------------------------------------------------------------------

void AHRS_VQF_SetTauAcc(AHRS_VQF_t* vqf, float tauAcc_s)
{
    if (vqf == NULL) {
        return;
    }

    if (vqf->params.tauAcc == tauAcc_s) {
        return;
    }
    vqf->params.tauAcc = tauAcc_s;

    double newB[3];
    double newA[2];

    vqf_filter_coeffs(vqf->params.tauAcc, vqf->coeffs.accTs, newB, newA);

    vqf_filter_adapt_state_for_coeff_change(
        vqf->state.lastAccLp, 3,
        vqf->coeffs.accLpB, vqf->coeffs.accLpA,
        newB, newA,
        vqf->state.accLpState);

    // For R and biasLP, the last output is not stored, approximate with state[0] (original approach).
    float R[9];
    for (size_t i = 0; i < 9; i++) {
        R[i] = (float)vqf->state.motionBiasEstRLpState[2 * i];
    }
    vqf_filter_adapt_state_for_coeff_change(
        R, 9,
        vqf->coeffs.accLpB, vqf->coeffs.accLpA,
        newB, newA,
        vqf->state.motionBiasEstRLpState);

    float biasLp[2];
    for (size_t i = 0; i < 2; i++) {
        biasLp[i] = (float)vqf->state.motionBiasEstBiasLpState[2 * i];
    }
    vqf_filter_adapt_state_for_coeff_change(
        biasLp, 2,
        vqf->coeffs.accLpB, vqf->coeffs.accLpA,
        newB, newA,
        vqf->state.motionBiasEstBiasLpState);

    memcpy(vqf->coeffs.accLpB, newB, sizeof(newB));
    memcpy(vqf->coeffs.accLpA, newA, sizeof(newA));
}

void AHRS_VQF_SetTauMag(AHRS_VQF_t* vqf, float tauMag_s)
{
    if (vqf == NULL) {
        return;
    }
    vqf->params.tauMag = tauMag_s;
    vqf->coeffs.kMag = vqf_gain_from_tau(vqf->params.tauMag, vqf->coeffs.magTs);
}

void AHRS_VQF_SetMotionBiasEstEnabled(AHRS_VQF_t* vqf, uint8_t enabled)
{
    if (vqf == NULL) {
        return;
    }
    if (vqf->params.motionBiasEstEnabled == enabled) {
        return;
    }
    vqf->params.motionBiasEstEnabled = enabled;
    vqf_fill_f64(vqf->state.motionBiasEstRLpState, 9 * 2, (double)VQF_NAN);
    vqf_fill_f64(vqf->state.motionBiasEstBiasLpState, 2 * 2, (double)VQF_NAN);
}

void AHRS_VQF_SetRestBiasEstEnabled(AHRS_VQF_t* vqf, uint8_t enabled)
{
    if (vqf == NULL) {
        return;
    }
    if (vqf->params.restBiasEstEnabled == enabled) {
        return;
    }
    vqf->params.restBiasEstEnabled = enabled;
    vqf->state.restDetected = 0;

    vqf_fill_f32(vqf->state.restLastSquaredDeviations, 2, 0.0f);
    vqf->state.restT = 0.0f;

    vqf_fill_f32(vqf->state.restLastGyrLp, 3, 0.0f);
    vqf_fill_f64(vqf->state.restGyrLpState, 3 * 2, (double)VQF_NAN);

    vqf_fill_f32(vqf->state.restLastAccLp, 3, 0.0f);
    vqf_fill_f64(vqf->state.restAccLpState, 3 * 2, (double)VQF_NAN);
}

void AHRS_VQF_SetMagDistRejectionEnabled(AHRS_VQF_t* vqf, uint8_t enabled)
{
    if (vqf == NULL) {
        return;
    }
    if (vqf->params.magDistRejectionEnabled == enabled) {
        return;
    }
    vqf->params.magDistRejectionEnabled = enabled;

    vqf->state.magDistDetected = 1;
    vqf->state.magRefNorm = 0.0f;
    vqf->state.magRefDip = 0.0f;
    vqf->state.magUndisturbedT = 0.0f;
    vqf->state.magRejectT = vqf->params.magMaxRejectionTime;
    vqf->state.magCandidateNorm = -1.0f;
    vqf->state.magCandidateDip = 0.0f;
    vqf->state.magCandidateT = 0.0f;
    vqf_fill_f64(vqf->state.magNormDipLpState, 2 * 2, (double)VQF_NAN);
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

void AHRS_VQF_SetBiasEstimate(AHRS_VQF_t* vqf, const float bias[3], float sigma)
{
    if (vqf == NULL) {
        return;
    }

    if (bias != NULL) {
        // Bias is provided in the external body frame convention
        vqf_vec_body_to_enu(bias, vqf->state.bias);
    }

    if (sigma > 0.0f) {
        const float P = vqf_square(sigma * (180.0f * 100.0f / VQF_PI_F));
        vqf_mat3_set_scaled_identity(P, vqf->state.biasP);
    }
}

// -----------------------------------------------------------------------------
// Diagnostics getters
// -----------------------------------------------------------------------------

uint8_t AHRS_VQF_GetRestDetected(const AHRS_VQF_t* vqf)
{
    return (vqf != NULL) ? vqf->state.restDetected : 0;
}

uint8_t AHRS_VQF_GetMagDistDetected(const AHRS_VQF_t* vqf)
{
    return (vqf != NULL) ? vqf->state.magDistDetected : 0;
}

void AHRS_VQF_GetRelativeRestDeviations(const AHRS_VQF_t* vqf, float out[2])
{
    if (vqf == NULL || out == NULL) {
        return;
    }

    out[0] = vqf_sqrt(vqf->state.restLastSquaredDeviations[0]) / (vqf->params.restThGyr * (VQF_PI_F / 180.0f));
    out[1] = vqf_sqrt(vqf->state.restLastSquaredDeviations[1]) / vqf->params.restThAcc;
}

float AHRS_VQF_GetMagRefNorm(const AHRS_VQF_t* vqf)
{
    return (vqf != NULL) ? vqf->state.magRefNorm : 0.0f;
}

float AHRS_VQF_GetMagRefDip(const AHRS_VQF_t* vqf)
{
    return (vqf != NULL) ? vqf->state.magRefDip : 0.0f;
}

// -----------------------------------------------------------------------------
// Bias estimate getter (ported from original, but with fixed includes / sizing)
// -----------------------------------------------------------------------------

float AHRS_VQF_GetBiasEstimate(const AHRS_VQF_t* vqf, float bias_out[3])
{
    if (vqf == NULL) {
        return 0.0f;
    }

    if (bias_out != NULL) {
        // Return bias in the external body frame convention
        vqf_vec_enu_to_body(vqf->state.bias, bias_out);
    }

    // Upper bound estimate for largest eigenvalue (Gershgorin circle theorem)
    const float sum1 = vqf_fabs(vqf->state.biasP[0]) + vqf_fabs(vqf->state.biasP[1]) + vqf_fabs(vqf->state.biasP[2]);
    const float sum2 = vqf_fabs(vqf->state.biasP[3]) + vqf_fabs(vqf->state.biasP[4]) + vqf_fabs(vqf->state.biasP[5]);
    const float sum3 = vqf_fabs(vqf->state.biasP[6]) + vqf_fabs(vqf->state.biasP[7]) + vqf_fabs(vqf->state.biasP[8]);
    const float lambda = vqf_max(sum1, vqf_max(sum2, sum3));

    // Convert from internal scaling back to rad/s and clip by biasSigmaInit
    const float sigma = vqf_sqrt(lambda) * (VQF_PI_F / (180.0f * 100.0f));
    return vqf_min(sigma, vqf->params.biasSigmaInit * (VQF_PI_F / 180.0f));
}

// -----------------------------------------------------------------------------
// Main update steps (ported from original VQF)
// -----------------------------------------------------------------------------

void AHRS_VQF_UpdateGyr(AHRS_VQF_t* vqf, const float gyr[3])
{
    if (vqf == NULL || gyr == NULL) {
        return;
    }

    // Convert input sample to internal ENU convention (VQF core runs in ENU)
    float gyr_enu[3];
    vqf_vec_body_to_enu(gyr, gyr_enu);
    gyr = gyr_enu;

    AHRS_VQF_Params_t* p = &vqf->params;
    AHRS_VQF_Coeffs_t* c = &vqf->coeffs;
    AHRS_VQF_State_t*  s = &vqf->state;

    // Rest detection (gyro)
    if (p->restBiasEstEnabled || p->magDistRejectionEnabled) {
        float tmp[3];
        vqf_filter_vec(gyr, 3, p->restFilterTau, c->gyrTs, c->restGyrLpB, c->restGyrLpA, s->restGyrLpState, tmp);
        memcpy(s->restLastGyrLp, tmp, sizeof(tmp));

        s->restLastSquaredDeviations[0] =
            vqf_square(gyr[0] - s->restLastGyrLp[0]) +
            vqf_square(gyr[1] - s->restLastGyrLp[1]) +
            vqf_square(gyr[2] - s->restLastGyrLp[2]);

        const float biasClip = p->biasClip * (VQF_PI_F / 180.0f);
        if (s->restLastSquaredDeviations[0] >= vqf_square(p->restThGyr * (VQF_PI_F / 180.0f)) ||
            vqf_fabs(s->restLastGyrLp[0]) > biasClip ||
            vqf_fabs(s->restLastGyrLp[1]) > biasClip ||
            vqf_fabs(s->restLastGyrLp[2]) > biasClip) {
            s->restT = 0.0f;
            s->restDetected = 0;
        }
    }

    // Remove estimated gyro bias
    const float gyrNoBias[3] = {gyr[0] - s->bias[0], gyr[1] - s->bias[1], gyr[2] - s->bias[2]};

    // Gyro prediction step
    const float gyrNorm = vqf_norm(gyrNoBias, 3);
    const float angle = gyrNorm * c->gyrTs;

    if (gyrNorm > VQF_EPS) {
        const float cA = vqf_cos(angle * 0.5f);
        const float sA = vqf_sin(angle * 0.5f) / gyrNorm;
        const float dq[4] = {cA, sA * gyrNoBias[0], sA * gyrNoBias[1], sA * gyrNoBias[2]};

        float qNew[4];
        vqf_quat_multiply(s->gyrQuat, dq, qNew);
        memcpy(s->gyrQuat, qNew, sizeof(qNew));
        vqf_normalize(s->gyrQuat, 4);
    }
}

void AHRS_VQF_UpdateAcc(AHRS_VQF_t* vqf, const float acc[3])
{
    if (vqf == NULL || acc == NULL) {
        return;
    }

    // Convert input sample to internal ENU convention (VQF core runs in ENU)
    float acc_enu[3];
    vqf_vec_body_to_enu(acc, acc_enu);
    acc = acc_enu;

    AHRS_VQF_Params_t* p = &vqf->params;
    AHRS_VQF_Coeffs_t* c = &vqf->coeffs;
    AHRS_VQF_State_t*  s = &vqf->state;

    // Ignore [0 0 0] samples
    if (acc[0] == 0.0f && acc[1] == 0.0f && acc[2] == 0.0f) {
        return;
    }

    // Rest detection (acc)
    if (p->restBiasEstEnabled) {
        float tmp[3];
        vqf_filter_vec(acc, 3, p->restFilterTau, c->accTs, c->restAccLpB, c->restAccLpA, s->restAccLpState, tmp);
        memcpy(s->restLastAccLp, tmp, sizeof(tmp));

        s->restLastSquaredDeviations[1] =
            vqf_square(acc[0] - s->restLastAccLp[0]) +
            vqf_square(acc[1] - s->restLastAccLp[1]) +
            vqf_square(acc[2] - s->restLastAccLp[2]);

        if (s->restLastSquaredDeviations[1] >= vqf_square(p->restThAcc)) {
            s->restT = 0.0f;
            s->restDetected = 0;
        } else {
            s->restT += c->accTs;
            if (s->restT >= p->restMinT) {
                s->restDetected = 1;
            }
        }
    }

    float accEarth[3];

    // Filter accelerometer in inertial frame
    vqf_quat_rotate(s->gyrQuat, acc, accEarth);
    vqf_filter_vec(accEarth, 3, p->tauAcc, c->accTs, c->accLpB, c->accLpA, s->accLpState, s->lastAccLp);

    // Transform to 6D earth frame and normalize
    vqf_quat_rotate(s->accQuat, s->lastAccLp, accEarth);
    vqf_normalize(accEarth, 3);

    // Inclination correction
    float accCorrQuat[4];
    const float q_w = vqf_sqrt((accEarth[2] + 1.0f) * 0.5f);
    if (q_w > 1e-6f) {
        accCorrQuat[0] = q_w;
        accCorrQuat[1] = 0.5f * accEarth[1] / q_w;
        accCorrQuat[2] = -0.5f * accEarth[0] / q_w;
        accCorrQuat[3] = 0.0f;
    } else {
        // Avoid numeric issues when acc is close to [0 0 -1]
        accCorrQuat[0] = 0.0f;
        accCorrQuat[1] = 1.0f;
        accCorrQuat[2] = 0.0f;
        accCorrQuat[3] = 0.0f;
    }

    float qNew[4];
    vqf_quat_multiply(accCorrQuat, s->accQuat, qNew);
    memcpy(s->accQuat, qNew, sizeof(qNew));
    vqf_normalize(s->accQuat, 4);

    // Correction angular rate (debug)
    s->lastAccCorrAngularRate = vqf_acos(accEarth[2]) / c->accTs;

    // Bias estimation
    if (p->motionBiasEstEnabled || p->restBiasEstEnabled) {
        const float biasClip = p->biasClip * (VQF_PI_F / 180.0f);

        float accGyrQuat[4];
        float R[9];
        float biasLp[2];

        // Rotation matrix from accGyrQuat (6D quaternion)
        vqf_get_quat6d_enu(vqf, accGyrQuat);

        R[0] = 1.0f - 2.0f * vqf_square(accGyrQuat[2]) - 2.0f * vqf_square(accGyrQuat[3]); // r11
        R[1] = 2.0f * (accGyrQuat[2] * accGyrQuat[1] - accGyrQuat[0] * accGyrQuat[3]);     // r12
        R[2] = 2.0f * (accGyrQuat[0] * accGyrQuat[2] + accGyrQuat[3] * accGyrQuat[1]);     // r13
        R[3] = 2.0f * (accGyrQuat[0] * accGyrQuat[3] + accGyrQuat[2] * accGyrQuat[1]);     // r21
        R[4] = 1.0f - 2.0f * vqf_square(accGyrQuat[1]) - 2.0f * vqf_square(accGyrQuat[3]); // r22
        R[5] = 2.0f * (accGyrQuat[2] * accGyrQuat[3] - accGyrQuat[1] * accGyrQuat[0]);     // r23
        R[6] = 2.0f * (accGyrQuat[3] * accGyrQuat[1] - accGyrQuat[0] * accGyrQuat[2]);     // r31
        R[7] = 2.0f * (accGyrQuat[0] * accGyrQuat[1] + accGyrQuat[3] * accGyrQuat[2]);     // r32
        R[8] = 1.0f - 2.0f * vqf_square(accGyrQuat[1]) - 2.0f * vqf_square(accGyrQuat[2]); // r33

        // Compute R*b_hat (x/y components)
        biasLp[0] = R[0] * s->bias[0] + R[1] * s->bias[1] + R[2] * s->bias[2];
        biasLp[1] = R[3] * s->bias[0] + R[4] * s->bias[1] + R[5] * s->bias[2];

        // Low-pass filter R and R*b_hat
        vqf_filter_vec(R, 9, p->tauAcc, c->accTs, c->accLpB, c->accLpA, s->motionBiasEstRLpState, R);
        vqf_filter_vec(biasLp, 2, p->tauAcc, c->accTs, c->accLpB, c->accLpA, s->motionBiasEstBiasLpState, biasLp);

        // Measurement error and covariance for Kalman update
        float w[3];
        float e[3];

        if (s->restDetected && p->restBiasEstEnabled) {
            e[0] = s->restLastGyrLp[0] - s->bias[0];
            e[1] = s->restLastGyrLp[1] - s->bias[1];
            e[2] = s->restLastGyrLp[2] - s->bias[2];
            vqf_mat3_set_scaled_identity(1.0f, R);
            w[0] = c->biasRestW;
            w[1] = c->biasRestW;
            w[2] = c->biasRestW;
        } else if (p->motionBiasEstEnabled) {
            e[0] = -accEarth[1] / c->accTs + biasLp[0] - (R[0] * s->bias[0] + R[1] * s->bias[1] + R[2] * s->bias[2]);
            e[1] =  accEarth[0] / c->accTs + biasLp[1] - (R[3] * s->bias[0] + R[4] * s->bias[1] + R[5] * s->bias[2]);
            e[2] = -(R[6] * s->bias[0] + R[7] * s->bias[1] + R[8] * s->bias[2]);
            w[0] = c->biasMotionW;
            w[1] = c->biasMotionW;
            w[2] = c->biasVerticalW;
        } else {
            w[0] = -1.0f; w[1] = -1.0f; w[2] = -1.0f;
        }

        // Step 1: P = P + V (increase covariance even if there is no measurement update)
        if (s->biasP[0] < c->biasP0) { s->biasP[0] += c->biasV; }
        if (s->biasP[4] < c->biasP0) { s->biasP[4] += c->biasV; }
        if (s->biasP[8] < c->biasP0) { s->biasP[8] += c->biasV; }

        if (w[0] >= 0.0f) {
            // Clip disagreement
            vqf_clip(e, 3, -biasClip, biasClip);

            // Step 2: K = P R^T inv(W + R P R^T)
            float K[9];

            vqf_mat3_multiply_tps_second(s->biasP, R, K);  // K = P R^T
            vqf_mat3_multiply(R, K, K);                    // K = R P R^T
            K[0] += w[0];
            K[4] += w[1];
            K[8] += w[2];                                  // K = W + R P R^T

            if (vqf_mat3_inv(K, K)) {                       // K = inv(...)
                vqf_mat3_multiply_tps_first(R, K, K);       // K = R^T inv(...)
                vqf_mat3_multiply(s->biasP, K, K);          // K = P R^T inv(...)

                // Step 3: bias = bias + K e
                s->bias[0] += K[0] * e[0] + K[1] * e[1] + K[2] * e[2];
                s->bias[1] += K[3] * e[0] + K[4] * e[1] + K[5] * e[2];
                s->bias[2] += K[6] * e[0] + K[7] * e[1] + K[8] * e[2];

                // Step 4: P = P - K R P
                vqf_mat3_multiply(K, R, K);                 // K = K R
                vqf_mat3_multiply(K, s->biasP, K);          // K = K R P
                for (size_t i = 0; i < 9; i++) {
                    s->biasP[i] -= K[i];
                }

                // Clip bias
                vqf_clip(s->bias, 3, -biasClip, biasClip);
            }
        }
    }
}

void AHRS_VQF_UpdateMag(AHRS_VQF_t* vqf, const float mag[3])
{
    if (vqf == NULL || mag == NULL) {
        return;
    }

    // Convert input sample to internal ENU convention (VQF core runs in ENU)
    float mag_enu[3];
    vqf_vec_body_to_enu(mag, mag_enu);
    mag = mag_enu;
// ignore [0 0 0] samples
    if (mag[0] == 0.0f && mag[1] == 0.0f && mag[2] == 0.0f) {
        return;
    }

    AHRS_VQF_Params_t* p = &vqf->params;
    AHRS_VQF_Coeffs_t* c = &vqf->coeffs;
    AHRS_VQF_State_t*  s = &vqf->state;

    float magEarth[3];

    // bring magnetometer measurement into 6D earth frame
    float accGyrQuat[4];
    vqf_get_quat6d_enu(vqf, accGyrQuat);
    vqf_quat_rotate(accGyrQuat, mag, magEarth);

    if (p->magDistRejectionEnabled) {
        s->magNormDip[0] = vqf_norm(magEarth, 3);
        // dip angle: asin can be replaced by 2*atan(x / sqrt(1-x*x))
        s->magNormDip[1] = -vqf_asin(magEarth[2] / s->magNormDip[0]);

        if (p->magCurrentTau > 0.0f) {
            vqf_filter_vec(s->magNormDip, 2, p->magCurrentTau, c->magTs,
                           c->magNormDipLpB, c->magNormDipLpA,
                           s->magNormDipLpState, s->magNormDip);
        }

        // magnetic disturbance detection
        if (vqf_fabs(s->magNormDip[0] - s->magRefNorm) < p->magNormTh * s->magRefNorm &&
            vqf_fabs(s->magNormDip[1] - s->magRefDip) < p->magDipTh * (VQF_PI_F / 180.0f)) {

            s->magUndisturbedT += c->magTs;
            if (s->magUndisturbedT >= p->magMinUndisturbedTime) {
                s->magDistDetected = 0;
                s->magRefNorm += c->kMagRef * (s->magNormDip[0] - s->magRefNorm);
                s->magRefDip  += c->kMagRef * (s->magNormDip[1] - s->magRefDip);
            }
        } else {
            s->magUndisturbedT = 0.0f;
            s->magDistDetected = 1;
        }

        // new magnetic field acceptance
        if (vqf_fabs(s->magNormDip[0] - s->magCandidateNorm) < p->magNormTh * s->magCandidateNorm &&
            vqf_fabs(s->magNormDip[1] - s->magCandidateDip) < p->magDipTh * (VQF_PI_F / 180.0f)) {

            if (vqf_norm(s->restLastGyrLp, 3) >= p->magNewMinGyr * (VQF_PI_F / 180.0f)) {
                s->magCandidateT += c->magTs;
            }

            s->magCandidateNorm += c->kMagRef * (s->magNormDip[0] - s->magCandidateNorm);
            s->magCandidateDip  += c->kMagRef * (s->magNormDip[1] - s->magCandidateDip);

            if (s->magDistDetected &&
                (s->magCandidateT >= p->magNewTime ||
                 (s->magRefNorm == 0.0f && s->magCandidateT >= p->magNewFirstTime))) {

                s->magRefNorm = s->magCandidateNorm;
                s->magRefDip  = s->magCandidateDip;
                s->magDistDetected = 0;
                s->magUndisturbedT = p->magMinUndisturbedTime;
            }
        } else {
            s->magCandidateT = 0.0f;
            s->magCandidateNorm = s->magNormDip[0];
            s->magCandidateDip  = s->magNormDip[1];
        }
    }

    // calculate disagreement angle based on current magnetometer measurement (ENU: atan2(E, N))
    s->lastMagDisAngle = vqf_atan2(magEarth[0], magEarth[1]) - s->delta;

    // make sure the disagreement angle is in the range [-pi, pi]
    if (s->lastMagDisAngle > VQF_PI_F) {
        s->lastMagDisAngle -= 2.0f * VQF_PI_F;
    } else if (s->lastMagDisAngle < -VQF_PI_F) {
        s->lastMagDisAngle += 2.0f * VQF_PI_F;
    }

    float k = c->kMag;

    if (p->magDistRejectionEnabled) {
        // magnetic disturbance rejection
        if (s->magDistDetected) {
            if (s->magRejectT <= p->magMaxRejectionTime) {
                s->magRejectT += c->magTs;
                k = 0.0f;
            } else {
                k /= p->magRejectionFactor;
            }
        } else {
            s->magRejectT = vqf_max(s->magRejectT - p->magRejectionFactor * c->magTs, 0.0f);
        }
    }

    // ensure fast initial convergence
    if (s->kMagInit != 0.0f) {
        // make sure that the gain k is at least 1/N, N=1,2,3,... in the first few samples
        if (k < s->kMagInit) {
            k = s->kMagInit;
        }

        // iterative expression to calculate 1/N
        s->kMagInit = s->kMagInit / (s->kMagInit + 1.0f);

        // disable if t > tauMag
        if (s->kMagInit * p->tauMag < c->magTs) {
            s->kMagInit = 0.0f;
        }
    }

    // first-order filter step
    s->delta += k * s->lastMagDisAngle;

    // calculate correction angular rate to facilitate debugging
    s->lastMagCorrAngularRate = k * s->lastMagDisAngle / c->magTs;

    // make sure delta is in the range [-pi, pi]
    if (s->delta > VQF_PI_F) {
        s->delta -= 2.0f * VQF_PI_F;
    } else if (s->delta < -VQF_PI_F) {
        s->delta += 2.0f * VQF_PI_F;
    }
}




