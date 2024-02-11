/* BEGIN Header */
/**
 ******************************************************************************
 * \file            IMU_Madgwick.c
 * \author          Andrea Vivani
 * \brief           Quaternion attitude estimator created by Sebastian Madgwick
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

/* Includes ------------------------------------------------------------------*/

#include "IMU_Madgwick.h"
#include "basicMath.h"
#include "math.h"
#include "quaternion.h"

/* Functions -----------------------------------------------------------------*/

void IMU_Madgwick_update(axis3f_t* angles, axis3f_t accel, axis3f_t gyro) {
    static quaternion_t q = {1, 0, 0, 0};
    //TODO make constants
    float beta = 0;
    beta = sqrtf(3.0f * 0.25f) * configIMU_MADGWICK_GYRO_ERROR;

    float inv_norm; //vector norm
    float SEqDot_omega_1, SEqDot_omega_2, SEqDot_omega_3,
        SEqDot_omega_4;                                       //quaternion derivative from gyroscope elements
    float f_1, f_2, f_3;                                      //objective function elements
    float J_11or24, J_12or23, J_13or22, J_14or21, J_32, J_33; //objective function Jacobian elements
    float SEqHatDot_1, SEqHatDot_2, SEqHatDot_3, SEqHatDot_4; //estimated direction of the gyroscope error

    /* Axulirary variables to avoid reapeated calculations */
    float halfSEq_1 = 0.5f * q.q0;
    float halfSEq_2 = 0.5f * q.q1;
    float halfSEq_3 = 0.5f * q.q2;
    float halfSEq_4 = 0.5f * q.q3;
    float twoSEq_1 = 2.0f * q.q0;
    float twoSEq_2 = 2.0f * q.q1;
    float twoSEq_3 = 2.0f * q.q2;

    /* Normalize the accelerometer measurement */
    inv_norm = -fastInvSqrt(accel.x * accel.x + accel.y * accel.y + accel.z * accel.z);
    if (isnan(inv_norm) || isinf(inv_norm)) {
        inv_norm = -1.f;
    }
    accel.x *= inv_norm;
    accel.y *= inv_norm;
    accel.z *= inv_norm;

    /* Compute the objective function and Jacobian */
    f_1 = twoSEq_2 * q.q3 - twoSEq_1 * q.q2 - accel.x;
    f_2 = twoSEq_1 * q.q1 + twoSEq_3 * q.q3 - accel.y;
    f_3 = 1.0f - twoSEq_2 * q.q1 - twoSEq_3 * q.q2 - accel.z;
    J_11or24 = twoSEq_3;
    J_12or23 = 2.0f * q.q3;
    J_13or22 = twoSEq_1;
    J_14or21 = twoSEq_2;
    J_32 = 2.0f * J_14or21;
    J_33 = 2.0f * J_11or24;

    /* Compute the gradient (matrix multiplication) */
    SEqHatDot_1 = J_14or21 * f_2 - J_11or24 * f_1;
    SEqHatDot_2 = J_12or23 * f_1 + J_13or22 * f_2 - J_32 * f_3;
    SEqHatDot_3 = J_12or23 * f_2 - J_33 * f_3 - J_13or22 * f_1;
    SEqHatDot_4 = J_14or21 * f_1 + J_11or24 * f_2;

    /* Normalise the gradient to estimate direction of the gyroscope error */
    inv_norm = fastInvSqrt(SEqHatDot_1 * SEqHatDot_1 + SEqHatDot_2 * SEqHatDot_2 + SEqHatDot_3 * SEqHatDot_3
                           + SEqHatDot_4 * SEqHatDot_4);
    if (isnan(inv_norm) || isinf(inv_norm)) {
        inv_norm = 1.f;
    }
    SEqHatDot_1 *= inv_norm;
    SEqHatDot_2 *= inv_norm;
    SEqHatDot_3 *= inv_norm;
    SEqHatDot_4 *= inv_norm;

    /* Compute the quaternion rate measured by gyroscopes */
    SEqDot_omega_1 = -halfSEq_2 * gyro.x - halfSEq_3 * gyro.y - halfSEq_4 * gyro.z;
    SEqDot_omega_2 = halfSEq_1 * gyro.x + halfSEq_3 * gyro.z - halfSEq_4 * gyro.y;
    SEqDot_omega_3 = halfSEq_1 * gyro.y - halfSEq_2 * gyro.z + halfSEq_4 * gyro.x;
    SEqDot_omega_4 = halfSEq_1 * gyro.z + halfSEq_2 * gyro.y - halfSEq_3 * gyro.x;

    /* Compute then integrate the estimated quaternion rate */
    q.q0 += (SEqDot_omega_1 - (beta * SEqHatDot_1)) * configIMU_MADGWICK_LOOP_TIME_S;
    q.q1 += (SEqDot_omega_2 - (beta * SEqHatDot_2)) * configIMU_MADGWICK_LOOP_TIME_S;
    q.q2 += (SEqDot_omega_3 - (beta * SEqHatDot_3)) * configIMU_MADGWICK_LOOP_TIME_S;
    q.q3 += (SEqDot_omega_4 - (beta * SEqHatDot_4)) * configIMU_MADGWICK_LOOP_TIME_S;

    /* Normalise quaternion */
    quaternionNorm(&q);

    /* Convert quaterionion to Euler angles */
    quaternionToEuler(&q, angles);

    return;
}
