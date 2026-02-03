/* BEGIN Header */
/**
 ******************************************************************************
 * \file            AHRS_EKF.h
 * \author          Andrea Vivani
 * \brief           Implementation of AV EKF for attitude and heading estimation
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
#ifndef __AHRS_EKF_H__
#define __AHRS_EKF_H__

#ifdef __cplusplus
extern "C" {
#endif
/* Includes ------------------------------------------------------------------*/

#include <stdint.h>
#include "commonTypes.h"
#include "main.h"

/* Configuration -------------------------------------------------------------*/

/* Initial value of the Roll angle, in rad */
#ifndef configAHRS_EKF_PHI0
#define configAHRS_EKF_PHI0 0.f
#endif

/* Initial value of the Pitch angle, in rad */
#ifndef configAHRS_EKF_THETA0
#define configAHRS_EKF_THETA0 0.f
#endif

/* Initial value of the Yaw angle, in rad */
#ifndef configAHRS_EKF_PSI0
#define configAHRS_EKF_PSI0 0.f
#endif

/* Initial value of the damping coefficient, inversely proportional to the mass of the quadrotor, in N*s/m */
#ifndef configAHRS_EKF_C_DAMP0
#define configAHRS_EKF_C_DAMP0 0.07413f
#endif

/* Initial value of the inclination of the magnetic field, in rad
Northern Emisphere positive (pointing down), Southern Emisphere negative (pointing up) */
#ifndef configAHRS_EKF_INCL0
#define configAHRS_EKF_INCL0 -1.0734f
#endif

/* Loop time, in s */
#ifndef configAHRS_EKF_LOOP_TIME_S
#error configAHRS_EKF_LOOP_TIME_S must be defined
#endif

/* Magnetometer correction loop time, in s */
#ifndef configAHRS_EKF_MAG_LOOP_TIME_S
#error configAHRS_EKF_MAG_LOOP_TIME_S must be defined
#endif

/* AHRS EKF noises */
/* Gyro x,y noise */
#ifndef configAHRS_EKF_GXY_NOISE
#define configAHRS_EKF_GXY_NOISE 1e-3f
#endif

/* Gyro z noise */
#ifndef configAHRS_EKF_GZ_NOISE
#define configAHRS_EKF_GZ_NOISE 1e-3f
#endif

/* Accel x,y noise */
#ifndef configAHRS_EKF_AXY_NOISE
#define configAHRS_EKF_AXY_NOISE 1e-2f
#endif

/* Accel z noise */
#ifndef configAHRS_EKF_AZ_NOISE
#define configAHRS_EKF_AZ_NOISE 1e-1f
#endif

/* Magnetometer noise */
#ifndef configAHRS_EKF_M_NOISE
#define configAHRS_EKF_M_NOISE 1e-2f
#endif

/* Damping coefficient noise */
#ifndef configAHRS_EKF_C_DAMP_NOISE
#define configAHRS_EKF_C_DAMP_NOISE 5e-4f
#endif

/* Magnetic field inclination noise */
#ifndef configAHRS_EKF_INCL_NOISE
#define configAHRS_EKF_INCL_NOISE 1e-3f
#endif

/* Bias acc z noise */
#ifndef configAHRS_EKF_B_AZ_NOISE
#define configAHRS_EKF_B_AZ_NOISE 1e-3f
#endif

/* Bias gyro noise */
#ifndef configAHRS_EKF_B_G_NOISE
#define configAHRS_EKF_B_G_NOISE 1e-4f
#endif

/* Velocity x,y local noise */
#ifndef configAHRS_EKF_VXY_NOISE
#define configAHRS_EKF_VXY_NOISE 1e-3f
#endif

/* Velocity z local noise */
#ifndef configAHRS_EKF_VZ_NOISE
#define configAHRS_EKF_VZ_NOISE 1e-3f
#endif

/* Velocity N,E global noise */
#ifndef configAHRS_EKF_VNE_NOISE
#define configAHRS_EKF_VNE_NOISE 1e-3f
#endif

/* Velocity d global noise (pointing down) */
#ifndef configAHRS_EKF_VD_NOISE
#define configAHRS_EKF_VD_NOISE 1e-3f
#endif

/* Function prototypes -------------------------------------------------------*/

/**
 * \brief           AHRS EKF filter initialization
 *
 * \param[out]      angles: Euler angles vector
 * \param[out]      velocities: translational velocities along local axes
 */
void AHRS_EKF_init(axis3f_t* angles, axis3f_t* velocities);

/**
 * \brief           Predict EKF state at current step
 *
 * \param[in]       az: accelerometer measurement along local z axis, in m/s^2
 * \param[in]       gyro: gyroscope measurements vector, in rad/s
 */
void AHRS_EKF_prediction(float az, axis3f_t gyro);

/**
 * \brief           Update EKF with accelerometer and gyro readings
 *
 * \param[out]      angles: Euler angles vector
 * \param[out]      velocities: translational velocities along local axes
 * \param[in]       accel: accelerometer measurements vector, in m/s^2
 * \param[in]       gyro: gyroscope measurements vector, in rad/s
 */
void AHRS_EKF_updateAccelGyro(axis3f_t* angles, axis3f_t* velocities, axis3f_t accel, axis3f_t gyro);

/**
 * \brief           Update EKF with magnetometer readings
 *
 * \param[out]      angles: Euler angles vector
 * \param[out]      velocities: translational velocities along local axes
 * \param[in]       mag: magnetometer measurements vector, in Gauss or mGauss
 */
void AHRS_EKF_updateMag(axis3f_t* angles, axis3f_t* velocities, axis3f_t mag);

/**
 * \brief           Update EKF with velocity readings along local x and y axis
 *
 * \param[out]      angles: Euler angles vector
 * \param[out]      velocities: translational velocities along local axes
 * \param[in]       vx: velocity measurement along local x axis, in m/s
 * \param[in]       vy: velocity measurement along local y axis, in m/s
 * \param[in]       dt_s: update loop time, in s
 */
void AHRS_EKF_updateVelXY(axis3f_t* angles, axis3f_t* velocities, float vx, float vy, float dt_s);

/**
 * \brief           Update EKF with velocity reading along local z axis
 *
 * \param[out]      angles: Euler angles vector
 * \param[out]      velocities: translational velocities along local axes
 * \param[in]       vz: velocity measurement along local z axis, in m/s
 * \param[in]       dt_s: update loop time, in s
 */
void AHRS_EKF_updateVelZ(axis3f_t* angles, axis3f_t* velocities, float vz, float dt_s);

/**
 * \brief           Update EKF with velocity readings along global N and E axis
 *
 * \param[out]      angles: Euler angles vector
 * \param[out]      velocities: translational velocities along local axes
 * \param[in]       vN: velocity measurement along global N axis, in m/s
 * \param[in]       vE: velocity measurement along global E axis, in m/s
 * \param[in]       dt_s: update loop time, in s
 */
void AHRS_EKF_updateVelNE(axis3f_t* angles, axis3f_t* velocities, float vN, float vE, float dt_s);

/**
 * \brief           Update EKF with velocity reading along global D axis
 *
 * \param[out]      angles: Euler angles vector
 * \param[out]      velocities: translational velocities along local axes
 * \param[in]       vD: velocity measurement along global D axis, in m/s
 * \param[in]       dt_s: update loop time, in s
 */
void AHRS_EKF_updateVelD(axis3f_t* angles, axis3f_t* velocities, float vD, float dt_s);

/**
 * \brief           Reset EKF to initial values
 * 
 * \param[out]      angles: Euler angles vector
 * \param[out]      velocities: translational velocities along local axes
 * \param[in]       phi0: initial roll value, in rad
 * \param[in]       theta0: initial pitch value, in rad
 * \param[in]       psi0: initial yaw value, in rad
 */
void AHRS_EKF_reset(axis3f_t* angles, axis3f_t* velocities, float phi0, float theta0, float psi0);

/**
 * \brief           Set EKF input noises
 * 
 * \param[in]       az: noise of z-axis acceleration
 * \param[in]       gxy: noise of x-y-axes rotational speed
 * \param[in]       gz: noise of z-axis rotational speed
 * \param[in]       c_damp: noise of translational damping coefficient
 * \param[in]       b_az: noise of z-axis acceleration bias
 * \param[in]       incl: noise of magnetic field inclination vector
 * \param[in]       b_g: noise of gyroscope bias
 */
void AHRS_EKF_setInputNoises(float gxy, float gz, float az, float c_damp, float b_az, float incl, float b_g);

/**
 * \brief           Set EKF measurement noises of acceleration along local x and y axis
 * 
 * \param[in]       axy: noise of x-y-axes acceleration measurement
 */
void AHRS_EKF_setAccelNoise(float axy);

/**
 * \brief           Set EKF measurement noises of magnetic field measurement
 * 
 * \param[in]       m: noise of magnetic field measurement
 */
void AHRS_EKF_setMagNoise(float m);

/**
 * \brief           Retrieve a specific value from the state
 *
 * \param[in]       idx: index of required element (Roll=Phi, Pitch=Theta, Yaw=Psi, Xd, Yd, Zd, c_damp, b_az, incl, b_gx, b_gy, b_gz)
 *
 * \return          Value of requested item
 */
float AHRS_EKF_getStateValue(uint8_t idx);

#ifdef __cplusplus
}
#endif

#endif /* __AHRS_EKF_H__ */