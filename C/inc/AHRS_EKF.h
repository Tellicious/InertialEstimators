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

/* Lower bound of the damping coefficient estimate, in N*s/m */
#ifndef configAHRS_EKF_C_DAMP_MIN
#define configAHRS_EKF_C_DAMP_MIN (0.1f * configAHRS_EKF_C_DAMP0)
#endif

/* Upper bound of the damping coefficient estimate, in N*s/m */
#ifndef configAHRS_EKF_C_DAMP_MAX
#define configAHRS_EKF_C_DAMP_MAX (20.f * configAHRS_EKF_C_DAMP0)
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
#ifndef configAHRS_EKF_MAG_UPDATE_TIME_S
#error configAHRS_EKF_MAG_UPDATE_TIME_S must be defined
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
#define configAHRS_EKF_C_DAMP_NOISE 1e-6f
#endif

/* Magnetic field inclination noise */
#ifndef configAHRS_EKF_INCL_NOISE
#define configAHRS_EKF_INCL_NOISE 1e-3f
#endif

/* Bias acc z noise */
#ifndef configAHRS_EKF_B_AZ_NOISE
#define configAHRS_EKF_B_AZ_NOISE 1e-4f
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

/* Initial state uncertainties (1-sigma), used to seed the covariance matrix on init and reset.
   These must reflect how wrong the initial state can actually be: seeding the covariance too small
   makes the filter overconfident and, since the attitude enters the model through sin() terms, it can
   settle on the inverted solution instead of the correct one. */

/* Initial standard deviation of the roll and pitch estimates, in rad */
#ifndef configAHRS_EKF_P0_ANGLES_STD
#define configAHRS_EKF_P0_ANGLES_STD 0.3f
#endif

/* Initial standard deviation of the yaw estimate, in rad. Kept separate from roll and pitch because yaw is
   only observable through the magnetometer, so its initial error is typically larger */
#ifndef configAHRS_EKF_P0_PSI_STD
#define configAHRS_EKF_P0_PSI_STD 0.5f
#endif

/* Initial standard deviation of the x,y local velocity estimates, in m/s */
#ifndef configAHRS_EKF_P0_VXY_STD
#define configAHRS_EKF_P0_VXY_STD 5.f
#endif

/* Initial standard deviation of the z local velocity estimate, in m/s */
#ifndef configAHRS_EKF_P0_VZ_STD
#define configAHRS_EKF_P0_VZ_STD 1.f
#endif

/* Initial standard deviation of the damping coefficient estimate, in N*s/m */
#ifndef configAHRS_EKF_P0_C_DAMP_STD
#define configAHRS_EKF_P0_C_DAMP_STD (0.5f * configAHRS_EKF_C_DAMP0)
#endif

/* Initial standard deviation of the z-axis acceleration bias estimate, in m/s^2 */
#ifndef configAHRS_EKF_P0_B_AZ_STD
#define configAHRS_EKF_P0_B_AZ_STD 0.5f
#endif

/* Initial standard deviation of the magnetic field inclination estimate, in rad */
#ifndef configAHRS_EKF_P0_INCL_STD
#define configAHRS_EKF_P0_INCL_STD 0.1f
#endif

/* Initial standard deviation of the gyroscope bias estimates, in rad/s */
#ifndef configAHRS_EKF_P0_B_G_STD
#define configAHRS_EKF_P0_B_G_STD 0.02f
#endif

/* Lower bound on the magnitude of cos(pitch), used to keep 1 / cos(pitch) finite near +-90 degrees */
#ifndef configAHRS_EKF_C_THETA_MIN
#define configAHRS_EKF_C_THETA_MIN 1e-3f
#endif

/* Number of elements of the EKF state vector */
#define AHRS_EKF_STATE_SIZE 12

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
 * \brief           Update EKF with accelerometer readings
 *
 * \param[out]      angles: Euler angles vector
 * \param[out]      velocities: translational velocities along local axes
 * \param[in]       accel: accelerometer measurements vector, in m/s^2
 */
void AHRS_EKF_updateAccel(axis3f_t* angles, axis3f_t* velocities, axis3f_t accel);

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
 * \param[in]       gxy: noise of x-y-axes rotational speed
 * \param[in]       gz: noise of z-axis rotational speed
 * \param[in]       az: noise of z-axis acceleration
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
 * \brief           Set EKF measurement noises of velocity along local x and y axis
 * 
 * \param[in]       vxy: noise of x-y-axes velocity measurement
 */
void AHRS_EKF_setVelXYNoise(float vxy);

/**
 * \brief           Set EKF measurement noises of velocity along local z axis
 * 
 * \param[in]       vz: noise of z-axis velocity measurement
 */
void AHRS_EKF_setVelZNoise(float vz);

/**
 * \brief           Set EKF measurement noises of velocity along global N and E axis
 * 
 * \param[in]       vNE: noise of N-E-axes velocity measurement
 */
void AHRS_EKF_setVelNENoise(float vNE);

/**
 * \brief           Set EKF measurement noises of velocity along global D axis
 * 
 * \param[in]       vD: noise of D-axis velocity measurement
 */
void AHRS_EKF_setVelDNoise(float vD); 

/**
 * \brief           Retrieve a specific value from the state
 *
 * \param[in]       idx: index of required element (Roll=Phi, Pitch=Theta, Yaw=Psi, Xd, Yd, Zd, c_damp, b_az, incl, b_gx, b_gy, b_gz)
 *
 * \return          Value of requested item, 0 if idx is out of range
 */
float AHRS_EKF_getStateValue(uint8_t idx);

#ifdef __cplusplus
}
#endif

#endif /* __AHRS_EKF_H__ */