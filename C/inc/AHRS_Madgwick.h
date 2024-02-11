/* BEGIN Header */
/**
 ******************************************************************************
 * \file            AHRS_Madgwick.h
 * \author          Andrea Vivani
 * \brief           Quaternion attitude and heading estimator created by
 *                  Sebastian Madgwick
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
#ifndef __AHRS_MADGWICK_H__
#define __AHRS_MADGWICK_H__

#ifdef __cplusplus
extern "C" {
#endif
/* Includes ------------------------------------------------------------------*/

#include <stdint.h>
#include "main.h"
#include "quaternion.h"

/* Configuration -------------------------------------------------------------*/
/* Gyroscope supposed measurement error (rad/s) */
#define configAHRS_MADGWICK_GYRO_ERROR 0.075574973510005f //sqrtf(3.0f * 0.25f) * 0.0872664626f
/* Gyroscope supposed bias drift (rad/s/s) */
#define configAHRS_MADGWICK_GYRO_DRIFT 0.003576860579028f //sqrtf(3.0f * 0.35f) * 0.0034906585f

#ifndef configAHRS_MADGWICK_LOOP_TIME_S
#error configAHRS_MADGWICK_LOOP_TIME_S must be defined
#endif

/* Function prototypes -------------------------------------------------------*/

/**
 * \brief           Update attitude estimation
 *
 * \param[out]      angles: Euler angles output vector
 * \param[in]       accel: accelerometer measurements vector, in m/s^2
 * \param[in]       gyro: gyroscope measurements vector, in rad/s
 * \param[in]       mag: magnetometer measurements vector, in Gauss or mGauss
 */
void AHRS_Madgwick_update(axis3f_t* angles, axis3f_t accel, axis3f_t gyro, axis3f_t mag);

#ifdef __cplusplus
}
#endif

#endif /* __AHRS_MADGWICK_H__ */