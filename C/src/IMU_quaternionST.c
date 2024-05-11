/* BEGIN Header */
/**
 ******************************************************************************
 * \file            IMU_quaternionST.c
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

/* Includes ------------------------------------------------------------------*/

#include "IMU_quaternionST.h"
#include <math.h>
#include "basicMath.h"
#include "quaternion.h"

/* Macros --------------------------------------------------------------------*/

#ifdef USE_FAST_MATH
#define SIN(x)     fastSin(x)
#define COS(x)     fastCos(x)
#define SQRT(x)    fastSqrt(x);
#define INVSQRT(x) fastInvSqrt(x);
#else
#define SIN(x)     sinf(x)
#define COS(x)     cosf(x)
#define SQRT(x)    sqrtf(x);
#define INVSQRT(x) 1.0f / sqrtf(x);
#endif /* USE_FAST_MATH */

/* Functions -----------------------------------------------------------------*/

void IMU_quaternionST_update(axis3f_t* angles, axis3f_t accel, axis3f_t gyro, IMU_KP_val_t kpVal) {
    const float halfT = 0.5f * configIMU_QUATERNIONST_LOOP_TIME_S;
    const float kInt = configIMU_QUATERNIONST_KI * configIMU_QUATERNIONST_LOOP_TIME_S;
    static quaternion_t q = {1, 0, 0, 0};
    static float exInt = 0, eyInt = 0, ezInt = 0;
    float ahrs_kp;
    float axf, ayf, azf, gxf, gyf, gzf;
    float norm;
    float vx, vy, vz;
    float ex, ey, ez;

    if (kpVal == IMU_KP_HIGH) {
        ahrs_kp = configIMU_QUATERNIONST_KP_HIGH;
    } else {
        ahrs_kp = configIMU_QUATERNIONST_KP_NORM;
    }

    axf = accel.x;
    ayf = accel.y;
    azf = accel.z;

#ifdef configIMU_QUATERNIONST_CORRECT_ACCEL_OFFSET
    axf += ((gyro.y * gyro.y + gyro.z * gyro.z) * configIMU_QUATERNIONST_ACCEL_OFFSET_X)
           - (gyro.x * gyro.y * configIMU_QUATERNIONST_ACCEL_OFFSET_Y)
           - (gyro.x * gyro.z * configIMU_QUATERNIONST_ACCEL_OFFSET_Z);
    ayf += ((gyro.x * gyro.x + gyro.z * gyro.z) * configIMU_QUATERNIONST_ACCEL_OFFSET_Y)
           - (gyro.x * gyro.y * configIMU_QUATERNIONST_ACCEL_OFFSET_X)
           - (gyro.y * gyro.z * configIMU_QUATERNIONST_ACCEL_OFFSET_Z);
    azf += ((gyro.x * gyro.x + gyro.y * gyro.y) * configIMU_QUATERNIONST_ACCEL_OFFSET_Z)
           - (gyro.x * gyro.y * configIMU_QUATERNIONST_ACCEL_OFFSET_X)
           - (gyro.y * gyro.z * configIMU_QUATERNIONST_ACCEL_OFFSET_Y);
#endif

    gxf = gyro.x;
    gyf = gyro.y;
    gzf = gyro.z;

    /* Normalise the accelerometer measurement */
    norm = INVSQRT(axf * axf + ayf * ayf + azf * azf);

    axf *= norm;
    ayf *= norm;
    azf *= norm;

    /* Estimate direction of gravity and flux (v and w) */
    /* NED reference frame. Accel reading when flat is [0; 0; -1] */
    vx = 2 * ((q.q0 * q.q2) - (q.q1 * q.q3));
    vy = -2 * ((q.q0 * q.q1) + (q.q2 * q.q3));
    // vz = ((q.q1 * q.q1) + (q.q2 * q.q2) - (q.q0 * q.q0) - (q.q3 * q.q3));
    vz = 2 * ((q.q1 * q.q1) + (q.q2 * q.q2)) - 1;

    ex = (ayf * vz - azf * vy);
    ey = (azf * vx - axf * vz);
    ez = (axf * vy - ayf * vx);

    /* Integral error scaled integral gain */
    exInt += ex * kInt;
    eyInt += ey * kInt;
    ezInt += ez * kInt;

    /* Adjust gyroscope measurements */
    gxf += ahrs_kp * ex + exInt;
    gyf += ahrs_kp * ey + eyInt;
    gzf += ahrs_kp * ez + ezInt;

    /* Integrate quaternion rate and normalise */
    q.q0 += (-q.q1 * gxf - q.q2 * gyf - q.q3 * gzf) * halfT;
    q.q1 += (q.q0 * gxf + q.q2 * gzf - q.q3 * gyf) * halfT;
    q.q2 += (q.q0 * gyf - q.q1 * gzf + q.q3 * gxf) * halfT;
    q.q3 += (q.q0 * gzf + q.q1 * gyf - q.q2 * gxf) * halfT;

    /* Normalise quaternion */
    quaternionNorm(&q);

    /* Convert quaterionion to Euler angles */
    quaternionToEuler(&q, angles);

    return;
}
