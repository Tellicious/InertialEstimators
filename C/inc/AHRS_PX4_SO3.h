/* BEGIN Header */
/**
 ******************************************************************************
 * \file            AHRS_PX4_SO3.h
 * \author          Andrea Vivani
 * \brief           Quaternion-based complementary filter based on Robert
 *                  Mahony work
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
#ifndef __AHRS_PX4_SO3_H__
#define __AHRS_PX4_SO3_H__

#ifdef __cplusplus
extern "C" {
#endif
/* Includes ------------------------------------------------------------------*/

#include <stdint.h>
#include "commonTypes.h"
#include "main.h"

/* Configuration -------------------------------------------------------------*/
/* Default Kp value */
#define configAHRS_PX4_SO3_KP 1
/* Default Ki value */
#define configAHRS_PX4_SO3_KI 0.05f

#ifndef configAHRS_PX4_SO3_LOOP_TIME_S
#error configAHRS_PX4_SO3_LOOP_TIME_S must be defined
#endif

/* Typedefs ------------------------------------------------------------------*/

/* Function prototypes -------------------------------------------------------*/

/**
 * \brief           Update attitude estimation
 *
 * \param[out]      angles: Euler angles output vector
 * \param[in]       accel: accelerometer measurements vector, in m/s^2
 * \param[in]       gyro: gyroscope measurements vector, in rad/s
 * \param[in]       mag: magnetometer measurements vector, in Gauss or mGauss
 */
void AHRS_PX4_S03_update(axis3f_t* angles, axis3f_t accel, axis3f_t gyro, axis3f_t mag);

/**
 * \brief           Reset filter to initial values
 * 
 * \param[out]      angles: Euler angles vector
 * \param[in]       phi0: initial roll value, in rad
 * \param[in]       theta0: initial pitch value, in rad
 * \param[in]       psi0: initial yaw value, in rad
 */
void AHRS_PX4_S03_reset(axis3f_t* angles, float phi0, float theta0, float psi0);

#ifdef __cplusplus
}
#endif

#endif /* __AHRS_PX4_SO3_H__ */
