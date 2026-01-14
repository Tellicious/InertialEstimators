/* BEGIN Header */
/**
 ******************************************************************************
 * \file            IMU_EKF.h
 * \author          Andrea Vivani
 * \brief           Implementation of AV EKF for attitude estimation (excl. Yaw)
 ******************************************************************************
 * \copyright
 *
 * Copyright 2022 Andrea Vivani
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
 * \attention AXES DIRECTIONS (SENSORS AND ESTIMATION): X POINTING FORWARD 
 * (ROLL, PHI), Y POINTING RIGHT (PITCH, THETA), Z POINTING DOWN (YAW, PSI)
 * IF SENSORS HAVE DIFFERENT ORIENTATION, ROTATE READINGS PRIOR THAN INPUT 
 * THEM TO THE ESTIMATOR
 */
/* END Header */

/* Define to prevent recursive inclusion -------------------------------------*/
#ifndef __IMU_EKF_H__
#define __IMU_EKF_H__

#ifdef __cplusplus
extern "C" {
#endif
/* Includes ------------------------------------------------------------------*/

#include <stdint.h>
#include "commonTypes.h"
#include "main.h"

/* Configuration -------------------------------------------------------------*/

/* Initial value of the Roll angle, in rad */
#ifndef configIMU_EKF_PHI0
#define configIMU_EKF_PHI0 0.f
#endif

/* Initial value of the Pitch angle, in rad */
#ifndef configIMU_EKF_THETA0
#define configIMU_EKF_THETA0 0.f
#endif

/* Initial value of the damping coefficient, inversely proportional to the mass of the quadrotor, in N*s/m */
#ifndef configIMU_EKF_C_DAMP0
#define configIMU_EKF_C_DAMP0 0.07413f
#endif

/* Loop time, in s */
#ifndef configIMU_EKF_LOOP_TIME_S
#error configIMU_EKF_LOOP_TIME_S must be defined
#endif

/* IMU EKF noises */
/* Gyro x,y noise */
#ifndef configIMU_EKF_GXY_NOISE
#define configIMU_EKF_GXY_NOISE 1e-3f
#endif

/* Gyro z noise */
#ifndef configIMU_EKF_GZ_NOISE
#define configIMU_EKF_GZ_NOISE 1e-3f
#endif

/* Accel x,y noise */
#ifndef configIMU_EKF_AXY_NOISE
#define configIMU_EKF_AXY_NOISE 1e-2f
#endif

/* Accel z noise */
#ifndef configIMU_EKF_AZ_NOISE
#define configIMU_EKF_AZ_NOISE 1e-1f
#endif

/* Damping coefficient noise */
#ifndef configIMU_EKF_C_DAMP_NOISE
#define configIMU_EKF_C_DAMP_NOISE 5e-4f
#endif

/* Bias acc z noise */
#ifndef configIMU_EKF_B_AZ_NOISE
#define configIMU_EKF_B_AZ_NOISE 1e-2f
#endif

/* Velocity x,y local noise */
#ifndef configIMU_EKF_VXY_NOISE
#define configIMU_EKF_VXY_NOISE 1e-3f
#endif

/* Velocity z local noise */
#ifndef configIMU_EKF_VZ_NOISE
#define configIMU_EKF_VZ_NOISE 1e-3f
#endif

/* Velocity d global noise (pointing down) */
#ifndef configIMU_EKF_VD_NOISE
#define configIMU_EKF_VD_NOISE 1e-1f
#endif

/* Function prototypes -------------------------------------------------------*/

/**
 * \brief           IMU EKF filter initialization
 *
 * \param[out]      angles: Euler angles vector
 * \param[out]      velocities: translational velocities along local axes
 */
void IMU_EKF_init(axis3f_t* angles, axis3f_t* velocities);

/**
 * \brief           Predict EKF state at current step
 *
 * \param[in]       az: accelerometer measurement along local z axis, in m/s^2
 * \param[in]       gyro: gyroscope measurements vector, in rad/s
 */
void IMU_EKF_prediction(float az, axis3f_t gyro);

/**
 * \brief           Update EKF with accelerometer and gyro readings
 *
 * \param[out]      angles: Euler angles vector
 * \param[out]      velocities: translational velocities along local axes
 * \param[in]       accel: accelerometer measurements vector, in m/s^2
 * \param[in]       gyro: gyroscope measurements vector, in rad/s
 */
void IMU_EKF_updateAccelGyro(axis3f_t* angles, axis3f_t* velocities, axis3f_t accel, axis3f_t gyro);

/**
 * \brief           Update EKF with velocity readings along local x and y axis
 *
 * \param[out]      angles: Euler angles vector
 * \param[out]      velocities: translational velocities along local axes
 * \param[in]       vx: velocity measurement along local x axis, in m/s
 * \param[in]       vy: velocity measurement along local y axis, in m/s
 * \param[in]       dt_s: update loop time, in s
 */
void IMU_EKF_updateVelXY(axis3f_t* angles, axis3f_t* velocities, float vx, float vy, float dt_s);

/**
 * \brief           Update EKF with velocity reading along local z axis
 *
 * \param[out]      angles: Euler angles vector
 * \param[out]      velocities: translational velocities along local axes
 * \param[in]       vz: velocity measurement along local z axis, in m/s
 * \param[in]       dt_s: update loop time, in s
 */
void IMU_EKF_updateVelZ(axis3f_t* angles, axis3f_t* velocities, float vz, float dt_s);

/**
 * \brief           Update EKF with velocity reading along global D axis
 *
 * \param[out]      angles: Euler angles vector
 * \param[out]      velocities: translational velocities along local axes
 * \param[in]       vD: velocity measurement along global D axis, in m/s
 * \param[in]       dt_s: update loop time, in s
 */
void IMU_EKF_updateVelD(axis3f_t* angles, axis3f_t* velocities, float vD, float dt_s);

/**
 * \brief           Reset EKF to initial values
 * 
 * \param[out]      angles: Euler angles vector
 * \param[out]      velocities: translational velocities along local axes
 * \param[in]       phi0: initial roll value, in rad
 * \param[in]       theta0: initial pitch value, in rad
 */
void IMU_EKF_reset(axis3f_t* angles, axis3f_t* velocities, float phi0, float theta0);

/**
 * \brief           Set EKF input noises
 * 
 * \param[in]       az: noise of z-axis acceleration
 * \param[in]       gxy: noise of x-y-axes rotational speed
 * \param[in]       gz: noise of z-axis rotational speed
 * \param[in]       c_damp: noise of translational damping coefficient
 * \param[in]       b_az: noise of z-axis acceleration bias
 */
void IMU_EKF_setInputNoises(float az, float gxy, float gz, float c_damp, float b_az);

/**
 * \brief           Set EKF measurement noises of acceleration along local x and y axis
 * 
 * \param[in]       axy: noise of x-y-axes acceleration measurement
 */
void IMU_EKF_setAccelNoise(float axy);

/**
 * \brief           Set EKF measurement noises of velocity along local x and y axis
 * 
 * \param[in]       vxy: noise of x-y-axes velocity measurement
 */
void IMU_EKF_setVelXYNoise(float vxy);

/**
 * \brief           Set EKF measurement noises of velocity along local z axis
 * 
 * \param[in]       vz: noise of z-axis velocity measurement
 */
void IMU_EKF_setVelZNoise(float vz);

/**
 * \brief           Set EKF measurement noises of velocity along global D axis
 * 
 * \param[in]       vD: noise of D-axis velocity measurement
 */
void IMU_EKF_setVelDNoise(float vD);

/**
 * \brief           Retrieve a specific value from the state
 *
 * \param[in]       idx: index of required element (Roll=Phi, Pitch=Theta, Xd, Yd, Zd, c_damp, b_az)
 *
 * \return          Value of requested item
 */
float IMU_EKF_getStateValue(uint8_t idx);

#ifdef __cplusplus
}
#endif

#endif /* __IMU_EKF_H__ */
