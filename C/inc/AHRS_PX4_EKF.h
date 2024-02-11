/* BEGIN Header */
/**
 ******************************************************************************
 * \file            AHRS_PX4_EKF.h
 * \author          Andrea Vivani
 * \brief           PX4 attitude and heading EKF estimator
 ******************************************************************************
 * \copyright
 *
 * Copyright 2024 Andrea Vivani
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
#ifndef __AHRS_PX4_EKF_H__
#define __AHRS_PX4_EKF_H__

#ifdef __cplusplus
extern "C" {
#endif
/* Includes ------------------------------------------------------------------*/

#include <stdint.h>
#include "commonTypes.h"
#include "main.h"

/* Configuration -------------------------------------------------------------*/

/* Rotational speed process noise */
#ifndef configAHRS_PX4_EKF_R_S_NOISE
#define configAHRS_PX4_EKF_R_S_NOISE 2e-2f
#endif

/* Rotational acceleration process noise */
#ifndef configAHRS_PX4_EKF_R_A_NOISE
#define configAHRS_PX4_EKF_R_A_NOISE 16.f
#endif

/* Acceleration process noise */
#ifndef configAHRS_PX4_EKF_A_NOISE
#define configAHRS_PX4_EKF_A_NOISE 1.8f
#endif

/* Magnetic field process noise */
#ifndef configAHRS_PX4_EKF_M_NOISE
#define configAHRS_PX4_EKF_M_NOISE 1.f
#endif

/* Gyroscope noise */
#ifndef configAHRS_PX4_EKF_GYRO_NOISE
#define configAHRS_PX4_EKF_GYRO_NOISE 4e-6f
#endif

/* Accelerometer noise */
#ifndef configAHRS_PX4_EKF_ACCEL_NOISE
#define configAHRS_PX4_EKF_ACCEL_NOISE 50.f
#endif

/* Magnetometer noise */
#ifndef configAHRS_PX4_EKF_MAG_NOISE
#define configAHRS_PX4_EKF_MAG_NOISE 0.5f
#endif

/* Initial value of the inclination of the magnetic field, in rad */
/* Northern Emisphere positive (pointing down), Southern Emisphere negative (pointing up) */
#ifndef configAHRS_PX4_EKF_INCL0
#define configAHRS_PX4_EKF_INCL0 -1.0734f
#endif

/* Define whether to use diagonal or complete inertia matrix*/
#define configAHRS_PX4_EKF_COMPLETE_INERTIA_MATRIX 0
#define configAHRS_PX4_EKF_DIAGONAL_INERTIA_MATRIX 1

#ifndef configAHRS_PX4_EKF_INERTIA_MATRIX
#define configAHRS_PX4_EKF_INERTIA_MATRIX 2
#endif

/* Prediction loop time, in s */
#ifndef configAHRS_PX4_EKF_PRED_LOOP_TIME_S
#error configAHRS_PX4_EKF_PRED_LOOP_TIME_S must be defined
#endif

/* Accelerometer update loop time, in s */
#ifndef configAHRS_PX4_EKF_ACC_LOOP_TIME_S
#error configAHRS_PX4_EKF_ACC_LOOP_TIME_S must be defined
#endif

/* Gyroscope update loop time, in s */
#ifndef configAHRS_PX4_EKF_GYRO_LOOP_TIME_S
#error configAHRS_PX4_EKF_GYRO_LOOP_TIME_S must be defined
#endif

/* Magnetometer update loop time, in s */
#ifndef configAHRS_PX4_EKF_MAG_LOOP_TIME_S
#error configAHRS_PX4_EKF_MAG_LOOP_TIME_S must be defined
#endif

/* Function prototypes -------------------------------------------------------*/

/**
 * \brief           AHRS PX4 EKF filter initialization
 */
void AHRS_PX4_EKF_init();

/**
 * \brief           Predict EKF state at current step
 */
void AHRS_PX4_EKF_prediction();

/**
 * \brief           Update EKF with gyro readings
 *
 * \param[in]       gyro: gyroscope measurements vector, in rad/s
 */
void AHRS_PX4_EKF_updateGyro(axis3f_t gyro);

/**
 * \brief           Update EKF with accelerometer readings
 *
 * \param[in]       accel: accelerometer measurements vector, in m/s^2
 */
void AHRS_PX4_EKF_updateAccel(axis3f_t accel);

/**
 * \brief           Update EKF with magnetometer readings
 *
 * \param[in]       mag: magnetometer measurements vector, in Gauss or mGauss
 */
void AHRS_PX4_EKF_updateMag(axis3f_t mag);

/**
 * \brief           Calculate Euler angles
 *
 * \param[out]      angles: Euler angles vector
 */
void AHRS_PX4_EKF_calculateAngles(axis3f_t* angles);

/**
 * \brief           Set magnetic field inclination
 *
 * \param[in]      incl_angle: magnetic field inclination, in rad
 */
void AHRS_Attitude_PX4_EKF_setInclination(float incl_angle);

/**
 * \brief           Set EKF input noises
 * 
 * \param[in]       rot_sp: noise of rotational speed
 * \param[in]       rot_acc: noise of rotational accelerations
 * \param[in]       acc: noise of acceleration estimation
 * \param[in]       mag: noise of magnetic field estimation
 */
void AHRS_Attitude_PX4_EKF_setProcessNoise(float rot_sp, float rot_acc, float acc, float mag);

/**
 * \brief           Set EKF measurement noises of gyroscope measurement
 * 
 * \param[in]       g: noise of gyroscope measurement
 */
void AHRS_Attitude_PX4_EKF_setGyroNoise(float g);

/**
 * \brief           Set EKF measurement noises of acceleration measurement
 * 
 * \param[in]       a: noise of acceleration measurement
 */
void AHRS_Attitude_PX4_EKF_setAccelNoise(float a);

/**
 * \brief           Set EKF measurement noises of magnetic field measurement
 * 
 * \param[in]       m: noise of magnetic field measurement
 */
void AHRS_Attitude_PX4_EKF_setMagNoise(float m);

#ifdef __cplusplus
}
#endif

#endif /* __AHRS_PX4_EKF_H__ */
