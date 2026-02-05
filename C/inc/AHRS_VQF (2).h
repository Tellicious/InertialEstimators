/* BEGIN Header */
/**
 ******************************************************************************
 * \file            AHRS_VQF.h
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

/* Define to prevent recursive inclusion -------------------------------------*/
#ifndef __AHRS_VQF_H__
#define __AHRS_VQF_H__

#ifdef __cplusplus
extern "C" {
#endif

/* Includes ------------------------------------------------------------------*/

#include <stdint.h>

#include "commonTypes.h"
#include "matrix.h"
#include "quaternion.h"

/* Configuration -------------------------------------------------------------*/

#if !defined(ADVUTILS_USE_DYNAMIC_ALLOCATION)
#error ADVUTILS_USE_DYNAMIC_ALLOCATION must be set to use AHRS_VQF (dynamic matrix allocation)
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

    float   kMagRef;
    matrix_t magNormDipLpB;
    matrix_t magNormDipLpA;
} AHRS_VQF_Coeffs_t;
typedef struct {
    quaternion_t gyrQuat;
    quaternion_t accQuat;
    float        delta;

    uint8_t restDetected;
    uint8_t magDistDetected;

    axis3f_t lastAccLp;
    matrix_t accLpState;
    float    lastAccCorrAngularRate;

    float    kMagInit;
    float    lastMagDisAngle;
    float    lastMagCorrAngularRate;

    axis3f_t bias;
    matrix_t biasP;

    matrix_t motionBiasEstRLpState;
    matrix_t motionBiasEstBiasLpState;

    matrix_t restLastSquaredDeviations;
    float    restT;
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
typedef struct {
    AHRS_VQF_Params_t params;
    AHRS_VQF_Coeffs_t coeffs;
    AHRS_VQF_State_t  state;

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
} AHRS_VQF_t;

/* Function prototypes -------------------------------------------------------*/

/**
 * \brief           Initialize a VQF instance.
 *
 * \param[in,out]   vqf: pointer to the filter instance.
 * \param[in]       gyrTs_s: gyroscope sampling time [s].
 * \param[in]       accTs_s: accelerometer sampling time [s] (if <=0 it defaults to gyrTs_s).
 * \param[in]       magTs_s: magnetometer sampling time [s] (if <=0 it defaults to gyrTs_s).
 */
void AHRS_VQF_Init(AHRS_VQF_t* vqf, float gyrTs_s, float accTs_s, float magTs_s);

/**
 * \brief           Deinitialize a VQF instance.
 *
 * Frees all dynamically allocated matrices belonging to the instance.
 *
 * \param[in,out]   vqf: pointer to the filter instance.
 */
void AHRS_VQF_Deinit(AHRS_VQF_t* vqf);

/**
 * \brief           Reset the VQF state (keeps current parameters and coefficients).
 *
 * \param[in,out]   vqf: pointer to the filter instance.
 */
void AHRS_VQF_Reset(AHRS_VQF_t* vqf);

/**
 * \brief           Update the 3D (gyro-only) part of the filter.
 *
 * Input gyroscope vector is expressed in the body-NED frame.
 * Units: rad/s.
 *
 * \param[in,out]   vqf: pointer to the filter instance.
 * \param[in]       gyr_rad_s: gyroscope measurement [rad/s] in body-NED.
 */
void AHRS_VQF_UpdateGyr(AHRS_VQF_t* vqf, axis3f_t gyr_rad_s);

/**
 * \brief           Update the 6D (acc) part of the filter.
 *
 * Input accelerometer vector is expressed in the body-NED frame.
 * The filter uses only the direction (norm used for thresholds).
 *
 * \param[in,out]   vqf: pointer to the filter instance.
 * \param[in]       acc: accelerometer measurement in body-NED.
 */
void AHRS_VQF_UpdateAcc(AHRS_VQF_t* vqf, axis3f_t acc);

/**
 * \brief           Update the 9D (mag) part of the filter.
 *
 * Input magnetometer vector is expressed in the body-NED frame.
 *
 * \param[in,out]   vqf: pointer to the filter instance.
 * \param[in]       mag: magnetometer measurement in body-NED.
 */
void AHRS_VQF_UpdateMag(AHRS_VQF_t* vqf, axis3f_t mag);

/**
 * \brief           Get the 3D gyro-only orientation.
 *
 * The returned quaternion represents the rotation from body-NED to NED.
 *
 * \param[in]       vqf: pointer to the filter instance.
 * \param[out]      q_out: output quaternion (scalar-first).
 */
void AHRS_VQF_GetQuat3D(const AHRS_VQF_t* vqf, quaternion_t* q_out);

/**
 * \brief           Get the 6D (gyro + acc) orientation.
 *
 * The returned quaternion represents the rotation from body-NED to NED.
 *
 * \param[in]       vqf: pointer to the filter instance.
 * \param[out]      q_out: output quaternion (scalar-first).
 */
void AHRS_VQF_GetQuat6D(const AHRS_VQF_t* vqf, quaternion_t* q_out);

/**
 * \brief           Get the 9D (gyro + acc + mag) orientation.
 *
 * The returned quaternion represents the rotation from body-NED to NED.
 *
 * \param[in]       vqf: pointer to the filter instance.
 * \param[out]      q_out: output quaternion (scalar-first).
 */
void AHRS_VQF_GetQuat9D(const AHRS_VQF_t* vqf, quaternion_t* q_out);

/**
 * \brief           Get the estimated yaw correction delta.
 *
 * Returned in NED convention, i.e. positive rotation about +Down.
 *
 * \param[in]       vqf: pointer to the filter instance.
 *
 * \return          Yaw correction delta [rad].
 */
float AHRS_VQF_GetDelta(const AHRS_VQF_t* vqf);

/**
 * \brief           Get current gyro bias estimate.
 *
 * Bias is expressed in the body-NED frame, same units as gyro input.
 *
 * \param[in]       vqf: pointer to the filter instance.
 * \param[out]      bias_out: if not NULL, receives the bias vector in body-NED.
 *
 * \return          A conservative bias sigma estimate [rad/s].
 */
float AHRS_VQF_GetBiasEstimate(const AHRS_VQF_t* vqf, axis3f_t* bias_out);

/**
 * \brief           Set the current gyro bias estimate.
 *
 * \param[in,out]   vqf: pointer to the filter instance.
 * \param[in]       bias: bias vector in body-NED [rad/s].
 * \param[in]       sigma: bias sigma [rad/s]. If <=0, covariance is left unchanged.
 */
void AHRS_VQF_SetBiasEstimate(AHRS_VQF_t* vqf, axis3f_t bias, float sigma);

/**
 * \brief           Return rest detection flag.
 *
 * \param[in]       vqf: pointer to the filter instance.
 *
 * \return          1 if rest detected, otherwise 0.
 */
uint8_t AHRS_VQF_GetRestDetected(const AHRS_VQF_t* vqf);

/**
 * \brief           Return magnetic disturbance detection flag.
 *
 * \param[in]       vqf: pointer to the filter instance.
 *
 * \return          1 if magnetic disturbance detected, otherwise 0.
 */
uint8_t AHRS_VQF_GetMagDistDetected(const AHRS_VQF_t* vqf);

/**
 * \brief           Set magnetic reference norm and dip.
 *
 * \param[in,out]   vqf: pointer to the filter instance.
 * \param[in]       norm: reference magnetic norm (same units as magnetometer).
 * \param[in]       dip_rad: reference dip angle [rad] (positive down).
 */
void AHRS_VQF_SetMagRef(AHRS_VQF_t* vqf, float norm, float dip_rad);

#ifdef __cplusplus
}
#endif

#endif /* __AHRS_VQF_H__ */
