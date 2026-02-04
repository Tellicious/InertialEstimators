/* BEGIN Header */
/**
 ******************************************************************************
 * \file            AHRS_VQF.h
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

#ifndef __AHRS_VQF_H__
#define __AHRS_VQF_H__

#ifdef __cplusplus
extern "C" {
#endif

/* Includes ------------------------------------------------------------------*/

#include <stdint.h>
#include "basicMath.h"
#include "matrix.h"
#include "quaternion.h"

/* Typedefs ------------------------------------------------------------------*/

/**
 * \brief VQF tuning parameters (see original VQF defaults).
 *
 * Notes:
 * - Time constants are expressed in seconds.
 * - Thresholds that represent angles or angular rates are expressed in degrees
 *   (as in the reference implementation) unless otherwise stated.
 */
typedef struct {
    /* Core filter time constants (seconds) */
    float tauAcc; /**< Accelerometer correction time constant [s] */
    float tauMag; /**< Magnetometer correction time constant [s] */

    /* Bias estimation / disturbance rejection enable flags */
    uint8_t motionBiasEstEnabled;    /**< Enable bias estimation during motion */
    uint8_t restBiasEstEnabled;      /**< Enable bias estimation during rest */
    uint8_t magDistRejectionEnabled; /**< Enable magnetic disturbance rejection */

    /* Gyro bias estimator tuning */
    float biasSigmaInit;                /**< Initial bias std-dev [deg/s] */
    float biasForgettingTime;           /**< Bias forgetting time [s] */
    float biasClip;                     /**< Bias clipping limit [deg/s] */
    float biasSigmaMotion;              /**< Motion bias std-dev [deg/s] */
    float biasVerticalForgettingFactor; /**< Vertical forgetting factor [-] */
    float biasSigmaRest;                /**< Rest bias std-dev [deg/s] */

    /* Rest detection */
    float restMinT;       /**< Minimum rest time [s] */
    float restFilterTau;  /**< Rest detection LP time constant [s] */
    float restThGyr;      /**< Gyro rest threshold [deg/s] */
    float restThAcc;      /**< Acc rest threshold [same as accel input] */

    /* Magnetic disturbance rejection */
    float magCurrentTau;          /**< Magnetic current LP time constant [s] */
    float magRefTau;              /**< Magnetic reference LP time constant [s] */
    float magNormTh;              /**< Magnetic norm threshold [-] */
    float magDipTh;               /**< Magnetic dip threshold [deg] */
    float magNewTime;             /**< New magnetic field accept time [s] */
    float magNewFirstTime;        /**< First mag accept time [s] */
    float magNewMinGyr;           /**< Minimum gyro rate for acceptance [deg/s] */
    float magMinUndisturbedTime;  /**< Minimum undisturbed time [s] */
    float magMaxRejectionTime;    /**< Maximum full rejection time [s] */
    float magRejectionFactor;     /**< Rejection factor [-] */
} AHRS_VQF_Params_t;

/**
 * \brief Internal coefficient cache (derived from parameters and sampling times).
 *
 * The coefficient vectors are stored as small column matrices for consistency with
 * the repository math utilities (matrix_t).
 */
typedef struct {
    float gyrTs; /**< Gyroscope sample period [s] */
    float accTs; /**< Accelerometer sample period [s] */
    float magTs; /**< Magnetometer sample period [s] */

    /* Accelerometer LP coefficients */
    matrix_t accLpB;        /**< 3x1 */
    float    accLpBData[3];
    matrix_t accLpA;        /**< 2x1 */
    float    accLpAData[2];

    float kMag; /**< Magnetometer correction gain */

    /* Bias estimator parameters */
    float biasP0;
    float biasV;
    float biasMotionW;
    float biasVerticalW;
    float biasRestW;

    /* Rest detection LP coefficients */
    matrix_t restGyrLpB;    /**< 3x1 */
    float    restGyrLpBData[3];
    matrix_t restGyrLpA;    /**< 2x1 */
    float    restGyrLpAData[2];

    matrix_t restAccLpB;    /**< 3x1 */
    float    restAccLpBData[3];
    matrix_t restAccLpA;    /**< 2x1 */
    float    restAccLpAData[2];

    /* Mag reference / current filter */
    float kMagRef;

    matrix_t magNormDipLpB; /**< 3x1 */
    float    magNormDipLpBData[3];
    matrix_t magNormDipLpA; /**< 2x1 */
    float    magNormDipLpAData[2];
} AHRS_VQF_Coeffs_t;

/**
 * \brief VQF runtime state.
 *
 * All internal computations run in ENU. Public APIs accept/return body-NED
 * quantities (axis3f_t) and return body->NED quaternions.
 */
typedef struct {
    quaternion_t gyrQuat; /**< 3D orientation from gyroscope integration (body->ENU) */
    quaternion_t accQuat; /**< Inclination correction quaternion (ENU) */
    float delta;          /**< Yaw correction angle about ENU +Z (Up) [rad] */

    uint8_t restDetected;
    uint8_t magDistDetected;

    axis3f_t lastAccLp;

    matrix_t accLpState;        /**< 6x1 (3 signals, 2 states each) */
    float    accLpStateData[3 * 2];

    float lastAccCorrAngularRate;

    float kMagInit;
    float lastMagDisAngle;
    float lastMagCorrAngularRate;

    axis3f_t bias; /**< Gyro bias estimate (body-ENU), same unit as gyro input */

    matrix_t biasP;             /**< 3x3 */
    float    biasPData[9];

    matrix_t motionBiasEstRLpState; /**< 18x1 (9 signals, 2 states each) */
    float    motionBiasEstRLpStateData[9 * 2];

    matrix_t motionBiasEstBiasLpState; /**< 4x1 (2 signals, 2 states each) */
    float    motionBiasEstBiasLpStateData[2 * 2];

    matrix_t restLastSquaredDeviations; /**< 2x1 */
    float    restLastSquaredDeviationsData[2];

    float restT;
    axis3f_t restLastGyrLp;

    matrix_t restGyrLpState;    /**< 6x1 */
    float    restGyrLpStateData[3 * 2];

    axis3f_t restLastAccLp;

    matrix_t restAccLpState;    /**< 6x1 */
    float    restAccLpStateData[3 * 2];

    float magRefNorm;
    float magRefDip;
    float magUndisturbedT;
    float magRejectT;

    float magCandidateNorm;
    float magCandidateDip;
    float magCandidateT;

    matrix_t magNormDip;        /**< 2x1 */
    float    magNormDipData[2];

    matrix_t magNormDipLpState; /**< 4x1 */
    float    magNormDipLpStateData[2 * 2];
} AHRS_VQF_State_t;

/**
 * \brief VQF filter instance.
 */
typedef struct {
    AHRS_VQF_Params_t params; /**< User parameters */
    AHRS_VQF_Coeffs_t coeffs; /**< Derived coefficients */
    AHRS_VQF_State_t  state;  /**< Runtime state */
} AHRS_VQF_t;

/* Function prototypes -------------------------------------------------------*/

/**
 * \brief Initialize a VQF instance.
 *
 * \param[in,out] vqf      VQF instance.
 * \param[in]     gyrTs_s  Gyroscope sample period [s].
 * \param[in]     accTs_s  Accelerometer sample period [s]. If <= 0, gyrTs_s is used.
 * \param[in]     magTs_s  Magnetometer sample period [s]. If <= 0, gyrTs_s is used.
 */
void AHRS_VQF_Init(AHRS_VQF_t* vqf, float gyrTs_s, float accTs_s, float magTs_s);

/**
 * \brief Reset the internal state of the filter (keeps current parameters and sampling times).
 *
 * \param[in,out] vqf VQF instance.
 */
void AHRS_VQF_Reset(AHRS_VQF_t* vqf);

/**
 * \brief Update the filter with a gyroscope sample.
 *
 * The input is assumed to be expressed in the body NED frame (x=N, y=E, z=D).
 *
 * \param[in,out] vqf       VQF instance.
 * \param[in]     gyr_rad_s Gyroscope sample [rad/s] (body-NED).
 */
void AHRS_VQF_UpdateGyr(AHRS_VQF_t* vqf, axis3f_t gyr_rad_s);

/**
 * \brief Update the filter with an accelerometer sample.
 *
 * The input is assumed to be expressed in the body NED frame (x=N, y=E, z=D).
 * Only the direction is used for inclination correction; the magnitude is used
 * only for rest detection thresholds.
 *
 * \param[in,out] vqf VQF instance.
 * \param[in]     acc Accelerometer sample (body-NED).
 */
void AHRS_VQF_UpdateAcc(AHRS_VQF_t* vqf, axis3f_t acc);

/**
 * \brief Update the filter with a magnetometer sample.
 *
 * The input is assumed to be expressed in the body NED frame (x=N, y=E, z=D).
 *
 * \param[in,out] vqf VQF instance.
 * \param[in]     mag Magnetometer sample (body-NED).
 */
void AHRS_VQF_UpdateMag(AHRS_VQF_t* vqf, axis3f_t mag);

/**
 * \brief Get the 3D orientation quaternion (gyro integration only).
 *
 * \param[in]  vqf   VQF instance.
 * \param[out] q_out Body->NED quaternion (scalar-first).
 */
void AHRS_VQF_GetQuat3D(const AHRS_VQF_t* vqf, quaternion_t* q_out);

/**
 * \brief Get the 6D orientation quaternion (gyro + accelerometer inclination correction).
 *
 * \param[in]  vqf   VQF instance.
 * \param[out] q_out Body->NED quaternion (scalar-first).
 */
void AHRS_VQF_GetQuat6D(const AHRS_VQF_t* vqf, quaternion_t* q_out);

/**
 * \brief Get the 9D orientation quaternion (gyro + acc + magnetometer yaw correction).
 *
 * \param[in]  vqf   VQF instance.
 * \param[out] q_out Body->NED quaternion (scalar-first).
 */
void AHRS_VQF_GetQuat9D(const AHRS_VQF_t* vqf, quaternion_t* q_out);

/**
 * \brief Get the yaw correction delta in NED convention.
 *
 * The internal delta is defined about ENU +Z (Up). This function returns the
 * equivalent correction about NED +Z (Down), i.e. the sign is inverted.
 *
 * \param[in] vqf VQF instance.
 * \return Yaw correction [rad] (positive about +Down).
 */
float AHRS_VQF_GetDelta(const AHRS_VQF_t* vqf);

/**
 * \brief Get the estimated gyroscope bias.
 *
 * \param[in]  vqf      VQF instance.
 * \param[out] bias_out Bias estimate [rad/s] in body-NED frame.
 * \return Conservative 1-sigma estimate [rad/s].
 */
float AHRS_VQF_GetBiasEstimate(const AHRS_VQF_t* vqf, axis3f_t* bias_out);

/**
 * \brief Return current rest detection status.
 *
 * \param[in] vqf VQF instance.
 * \return 1 if rest is detected, 0 otherwise.
 */
uint8_t AHRS_VQF_GetRestDetected(const AHRS_VQF_t* vqf);

/**
 * \brief Return current magnetic disturbance detection status.
 *
 * \param[in] vqf VQF instance.
 * \return 1 if a magnetic disturbance is detected, 0 otherwise.
 */
uint8_t AHRS_VQF_GetMagDistDetected(const AHRS_VQF_t* vqf);

/**
 * \brief Get normalized deviations used by the rest detector.
 *
 * Output is a 2x1 vector:
 * - out[0] = gyro deviation / restThGyr
 * - out[1] = acc  deviation / restThAcc
 *
 * \param[in]  vqf VQF instance.
 * \param[out] out 2x1 matrix to receive the deviations.
 */
void AHRS_VQF_GetRelativeRestDeviations(const AHRS_VQF_t* vqf, matrix_t* out);

/**
 * \brief Get current magnetic field reference norm.
 *
 * \param[in] vqf VQF instance.
 * \return Magnetic field norm reference.
 */
float AHRS_VQF_GetMagRefNorm(const AHRS_VQF_t* vqf);

/**
 * \brief Get current magnetic field reference dip angle (rad).
 *
 * \param[in] vqf VQF instance.
 * \return Magnetic field dip reference [rad].
 */
float AHRS_VQF_GetMagRefDip(const AHRS_VQF_t* vqf);

/**
 * \brief Set accelerometer correction time constant.
 *
 * \param[in,out] vqf       VQF instance.
 * \param[in]     tauAcc_s  Time constant [s].
 */
void AHRS_VQF_SetTauAcc(AHRS_VQF_t* vqf, float tauAcc_s);

/**
 * \brief Set magnetometer correction time constant.
 *
 * \param[in,out] vqf       VQF instance.
 * \param[in]     tauMag_s  Time constant [s].
 */
void AHRS_VQF_SetTauMag(AHRS_VQF_t* vqf, float tauMag_s);

/**
 * \brief Enable/disable motion bias estimation.
 *
 * \param[in,out] vqf     VQF instance.
 * \param[in]     enabled 0/1.
 */
void AHRS_VQF_SetMotionBiasEstEnabled(AHRS_VQF_t* vqf, uint8_t enabled);

/**
 * \brief Enable/disable rest bias estimation.
 *
 * \param[in,out] vqf     VQF instance.
 * \param[in]     enabled 0/1.
 */
void AHRS_VQF_SetRestBiasEstEnabled(AHRS_VQF_t* vqf, uint8_t enabled);

/**
 * \brief Enable/disable magnetic disturbance rejection.
 *
 * \param[in,out] vqf     VQF instance.
 * \param[in]     enabled 0/1.
 */
void AHRS_VQF_SetMagDistRejectionEnabled(AHRS_VQF_t* vqf, uint8_t enabled);

/**
 * \brief Set rest detection thresholds.
 *
 * \param[in,out] vqf   VQF instance.
 * \param[in]     thGyr Gyro threshold [deg/s].
 * \param[in]     thAcc Acc threshold [same as accel input].
 */
void AHRS_VQF_SetRestDetectionThresholds(AHRS_VQF_t* vqf, float thGyr, float thAcc);

/**
 * \brief Set the magnetic field reference.
 *
 * \param[in,out] vqf  VQF instance.
 * \param[in]     norm Magnetic field norm.
 * \param[in]     dip  Magnetic field dip [rad].
 */
void AHRS_VQF_SetMagRef(AHRS_VQF_t* vqf, float norm, float dip);

/**
 * \brief Set the gyroscope bias estimate and covariance.
 *
 * \param[in,out] vqf      VQF instance.
 * \param[in]     bias     Bias estimate [rad/s] in body-NED frame.
 * \param[in]     sigma    1-sigma value [rad/s]. If <= 0, covariance is not changed.
 */
void AHRS_VQF_SetBiasEstimate(AHRS_VQF_t* vqf, axis3f_t bias, float sigma);

#ifdef __cplusplus
}
#endif

#endif /* __AHRS_VQF_H__ */
