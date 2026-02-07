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
 * all copies or substantial portions of the Software
 *
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE
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
#include "main.h"
#include "matrix.h"
#include "quaternion.h"

/* Configuration -------------------------------------------------------------*/

#if !defined(ADVUTILS_USE_DYNAMIC_ALLOCATION)
#error ADVUTILS_USE_DYNAMIC_ALLOCATION must be set to use AHRS_VQF (dynamic matrix allocation)
#endif

/* Loop time, in s */
#ifndef configAHRS_VQF_LOOP_TIME_S
#error configAHRS_VQF_LOOP_TIME_S must be defined
#endif

/* Magnetometer correction loop time, in s */
#ifndef configAHRS_VQF_MAG_LOOP_TIME_S
#error configAHRS_VQF_MAG_LOOP_TIME_S must be defined
#endif

/* Function prototypes -------------------------------------------------------*/

/**
 * \brief           Initialize a VQF instance
 */
void AHRS_VQF_Init();

/**
 * \brief           Deinitialize a VQF instance
 *
 * Frees all dynamically allocated matrices belonging to the instance
 */
void AHRS_VQF_Deinit();

/**
 * \brief           Reset the VQF state (keeps current parameters and coefficients)
 */
void AHRS_VQF_Reset();

/**
 * \brief           Update the 3D (gyro-only) part of the filter
 *
 * Input gyroscope vector is expressed in the body-NED frame
 * Units: rad/s
 *
 * \param[in]       gyro: gyroscope measurement [rad/s] in body-NED
 */
void AHRS_VQF_updateGyro(axis3f_t gyro);

/**
 * \brief           Update the 6D (acc) part of the filter
 *
 * Input accelerometer vector is expressed in the body-NED frame
 * The filter uses only the direction (norm used for thresholds)
 *
 * \param[in]       acc: accelerometer measurement in body-NED
 */
void AHRS_VQF_updateAcc(axis3f_t acc);

/**
 * \brief           Update the 9D (mag) part of the filter
 *
 * Input magnetometer vector is expressed in the body-NED frame
 *
 * \param[in]       mag: magnetometer measurement in body-NED
 */
void AHRS_VQF_updateMag(axis3f_t mag);

/**
 * \brief           Get the 6D (gyro + acc) orientation
 *
 * The returned quaternion represents the rotation from body-NED to NED
 *
 * \param[out]      angles: output Euler angles
 */
void AHRS_VQF_Get6D(axis3f_t* angles);

/**
 * \brief           Get the 9D (gyro + acc + mag) orientation
 *
 * The returned quaternion represents the rotation from body-NED to NED
 *
 * \param[out]      angles: output Euler angles
 */
void AHRS_VQF_Get9D(axis3f_t* angles);

/**
 * \brief           Get the estimated yaw correction delta
 *
 * Returned in NED convention, i.e. positive rotation about +Down
 *
 *
 * \return          Yaw correction delta [rad]
 */
float AHRS_VQF_GetDelta();

/**
 * \brief           Get current gyro bias estimate
 *
 * Bias is expressed in the body-NED frame, same units as gyro input
 *
 * \param[out]      bias_out: if not NULL, receives the bias vector in body-NED
 *
 * \return          A conservative bias sigma estimate [rad/s]
 */
float AHRS_VQF_GetBiasEstimate(axis3f_t* bias_out);

/**
 * \brief           Set the current gyro bias estimate
 *
 * \param[in]       bias: bias vector in body-NED [rad/s]
 * \param[in]       sigma: bias sigma [rad/s]. If <=0, covariance is left unchanged
 */
void AHRS_VQF_SetBiasEstimate(axis3f_t bias, float sigma);

/**
 * \brief           Return rest detection flag
 *
 *
 * \return          1 if rest detected, otherwise 0
 */
uint8_t AHRS_VQF_GetRestDetected();

/**
 * \brief           Return magnetic disturbance detection flag
 *
 *
 * \return          1 if magnetic disturbance detected, otherwise 0
 */
uint8_t AHRS_VQF_GetMagDistDetected();

/**
 * \brief           Set magnetic reference norm and dip
 *
 * \param[in]       norm: reference magnetic norm (same units as magnetometer)
 * \param[in]       dip_rad: reference dip angle [rad] (positive down)
 */
void AHRS_VQF_SetMagRef(float norm, float dip_rad);

#ifdef __cplusplus
}
#endif

#endif /* __AHRS_VQF_H__ */
