/* BEGIN Header */
/**
 ******************************************************************************
 * \file            IMU_quaternionST.h
 * \author          Andrea Vivani
 * \brief           Quaternion-based IMU from ST Drone
 ******************************************************************************
 * \copyright
 *
 * Copyright 2023 Andrea Vivani
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
#ifndef __IMU_QUATERNIONST_H__
#define __IMU_QUATERNIONST_H__

#ifdef __cplusplus
extern "C" {
#endif
/* Includes ------------------------------------------------------------------*/

#include <stdint.h>
#include "commonTypes.h"
#include "main.h"

/* Configuration -------------------------------------------------------------*/

#ifndef configIMU_QUATERNIONST_KP_HIGH
#define configIMU_QUATERNIONST_KP_HIGH 10.0f
#endif

#ifndef configIMU_QUATERNIONST_KP_NORM
#define configIMU_QUATERNIONST_KP_NORM 0.4f
#endif

#ifndef configIMU_QUATERNIONST_KI
#define configIMU_QUATERNIONST_KI 0.1f
#endif

#ifndef configIMU_QUATERNIONST_LOOP_TIME_S
#error configIMU_QUATERNIONST_LOOP_TIME_S must be defined
#endif

/* Typedefs ------------------------------------------------------------------*/

typedef enum { IMU_KP_HIGH = 0, IMU_KP_NORM = 1 } IMU_KP_val_t;

/* Function prototypes -------------------------------------------------------*/

/**
 * \brief           Update attitude estimation
 *
 * \param[out]      angles: Euler angles output vector
 * \param[in]       accel: accelerometer measurements vector, in m/s^2
 * \param[in]       gyro: gyroscope measurements vector, in rad/s
 * \param[in]       kpVal: parameter to define wether to use high (IMU_KP_HIGH) or normal (IMU_KP_NORM) kp values
 */
void IMU_quaternionST_update(axis3f_t* angles, axis3f_t accel, axis3f_t gyro, IMU_KP_val_t kpVal);

#ifdef __cplusplus
}
#endif

#endif /* __IMU_QUATERNIONST_H__ */
