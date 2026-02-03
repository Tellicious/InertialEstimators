// Copyright (c) 2024 Hugo Chiang
// SPDX-License-Identifier: MIT
//
// VQF (Versatile Quaternion-based Filter) integration for InertialEstimators
//
// Coordinate convention
// - Input samples (gyro/acc/mag) are expressed in the body frame.
// - This module is configured to accept body-NED IMU axes by default:
//     X: North / Forward
//     Y: East  / Right
//     Z: Down
// - Output quaternions represent the rotation from body frame to NED navigation frame.
//   NED axes: X North, Y East, Z Down.
// - Quaternion format: [w, x, y, z] (scalar-first).

#ifndef AHRS_VQF_H
#define AHRS_VQF_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>

/* Configuration -------------------------------------------------------------*/

/* Define the body input axes convention.
 *  - 1: inputs are in body-NED (X North/Forward, Y East/Right, Z Down)
 *  - 0: inputs are in body-ENU (X East, Y North, Z Up)
 */
#ifndef configAHRS_VQF_INPUT_BODY_NED
#define configAHRS_VQF_INPUT_BODY_NED 1
#endif


typedef struct
{
    // --- Core filter time constants (seconds) ---
    float tauAcc;
    float tauMag;

    // --- Bias estimation / disturbance rejection ---
    uint8_t  motionBiasEstEnabled;
    uint8_t  restBiasEstEnabled;
    uint8_t  magDistRejectionEnabled;

    // --- Gyro bias estimator tuning (units mostly follow original VQF defaults) ---
    float biasSigmaInit;
    float biasForgettingTime;
    float biasClip;
    float biasSigmaMotion;
    float biasVerticalForgettingFactor;
    float biasSigmaRest;

    // --- Rest detection ---
    float restMinT;
    float restFilterTau;
    float restThGyr;
    float restThAcc;

    // --- Magnetic disturbance rejection ---
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

typedef struct
{
    float  gyrTs;
    float  accTs;
    float  magTs;

    double accLpB[3];
    double accLpA[2];
    float  kMag;

    float  biasP0;
    float  biasV;
    float  biasMotionW;
    float  biasVerticalW;
    float  biasRestW;

    double restGyrLpB[3];
    double restGyrLpA[2];
    double restAccLpB[3];
    double restAccLpA[2];

    float  kMagRef;
    double magNormDipLpB[3];
    double magNormDipLpA[2];
} AHRS_VQF_Coeffs_t;

typedef struct
{
    float gyrQuat[4];
    float accQuat[4];
    float delta;

    uint8_t  restDetected;
    uint8_t  magDistDetected;

    float  lastAccLp[3];
    double accLpState[3 * 2];
    float  lastAccCorrAngularRate;

    float  kMagInit;
    float  lastMagDisAngle;
    float  lastMagCorrAngularRate;

    float bias[3];
    float biasP[9];

    double motionBiasEstRLpState[9 * 2];
    double motionBiasEstBiasLpState[2 * 2];

    float restLastSquaredDeviations[2];
    float restT;
    float restLastGyrLp[3];
    double restGyrLpState[3 * 2];
    float restLastAccLp[3];
    double restAccLpState[3 * 2];

    float magRefNorm;
    float magRefDip;
    float magUndisturbedT;
    float magRejectT;

    float magCandidateNorm;
    float magCandidateDip;
    float magCandidateT;

    float magNormDip[2];
    double magNormDipLpState[2 * 2];
} AHRS_VQF_State_t;

typedef struct
{
    AHRS_VQF_Params_t params;
    AHRS_VQF_Coeffs_t coeffs;
    AHRS_VQF_State_t  state;
} AHRS_VQF_t;

// -----------------------------------------------------------------------------
// Lifecycle
// -----------------------------------------------------------------------------

void AHRS_VQF_Init(AHRS_VQF_t* vqf, float gyrTs_s, float accTs_s, float magTs_s);
void AHRS_VQF_Reset(AHRS_VQF_t* vqf);

// -----------------------------------------------------------------------------
// Update functions
//  - Call these at the corresponding sensor sampling rate.
//  - Gyro units: rad/s.
//  - Acc units: m/s^2 or g (only direction is used; magnitude is used for thresholds).
//  - Mag units: arbitrary, but consistent.
// -----------------------------------------------------------------------------

void AHRS_VQF_UpdateGyr(AHRS_VQF_t* vqf, const float gyr[3]);
void AHRS_VQF_UpdateAcc(AHRS_VQF_t* vqf, const float acc[3]);
void AHRS_VQF_UpdateMag(AHRS_VQF_t* vqf, const float mag[3]);

// -----------------------------------------------------------------------------
// Outputs
// -----------------------------------------------------------------------------

// Body -> NED quaternion (scalar-first [w, x, y, z])
void AHRS_VQF_GetQuat3D(const AHRS_VQF_t* vqf, float q_out[4]);
void AHRS_VQF_GetQuat6D(const AHRS_VQF_t* vqf, float q_out[4]);
void AHRS_VQF_GetQuat9D(const AHRS_VQF_t* vqf, float q_out[4]);

// Yaw correction delta (radians) in NED convention (positive rotation about +Down).
float AHRS_VQF_GetDelta(const AHRS_VQF_t* vqf);

// Estimated gyro bias (same units as gyro input). Returns a conservative sigma estimate (same unit).
float AHRS_VQF_GetBiasEstimate(const AHRS_VQF_t* vqf, float bias_out[3]);

// Convenience / diagnostics
uint8_t  AHRS_VQF_GetRestDetected(const AHRS_VQF_t* vqf);
uint8_t  AHRS_VQF_GetMagDistDetected(const AHRS_VQF_t* vqf);
void  AHRS_VQF_GetRelativeRestDeviations(const AHRS_VQF_t* vqf, float out[2]);
float AHRS_VQF_GetMagRefNorm(const AHRS_VQF_t* vqf);
float AHRS_VQF_GetMagRefDip(const AHRS_VQF_t* vqf);

// -----------------------------------------------------------------------------
// Parameter setters (optional)
// -----------------------------------------------------------------------------

void AHRS_VQF_SetTauAcc(AHRS_VQF_t* vqf, float tauAcc_s);
void AHRS_VQF_SetTauMag(AHRS_VQF_t* vqf, float tauMag_s);
void AHRS_VQF_SetMotionBiasEstEnabled(AHRS_VQF_t* vqf, uint8_t enabled);
void AHRS_VQF_SetRestBiasEstEnabled(AHRS_VQF_t* vqf, uint8_t enabled);
void AHRS_VQF_SetMagDistRejectionEnabled(AHRS_VQF_t* vqf, uint8_t enabled);
void AHRS_VQF_SetRestDetectionThresholds(AHRS_VQF_t* vqf, float thGyr, float thAcc);
void AHRS_VQF_SetMagRef(AHRS_VQF_t* vqf, float norm, float dip);
void AHRS_VQF_SetBiasEstimate(AHRS_VQF_t* vqf, const float bias[3], float sigma);

#ifdef __cplusplus
}
#endif

#endif // AHRS_VQF_H
