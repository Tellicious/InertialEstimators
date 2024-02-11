/* BEGIN Header */
/**
 ******************************************************************************
 * \file            AHRS_Madgwick.c
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

/* Includes ------------------------------------------------------------------*/

#include "AHRS_Madgwick.h"
#include "basicMath.h"
#include "math.h"

/* Functions -----------------------------------------------------------------*/
void AHRS_Madgwick_update(axis3f_t* angles, axis3f_t accel, axis3f_t gyro, axis3f_t mag) { // local system variables
    static quaternion_t q = {1, 0, 0, 0};
    /* Initial magnetic field estimated direction */
    static float b_x = 1;
    static float b_z = 0;
    /* Initial gyroscope estimated biases */
    static axis3f_t w_b = {0, 0, 0};
    /* Filter parameters */
    float inv_norm;                                                       //inverse of vector norm
    float SEqDot_omega_1, SEqDot_omega_2, SEqDot_omega_3, SEqDot_omega_4; //quaternion rate from gyroscopes elements
    float f_1, f_2, f_3, f_4, f_5, f_6;                                   //objective function elements
    float J_11or24, J_12or23, J_13or22, J_14or21, J_32, J_33,             //objective function Jacobian elements
        J_41, J_42, J_43, J_44, J_51, J_52, J_53, J_54, J_61, J_62, J_63, J_64;
    float SEqHatDot_1, SEqHatDot_2, SEqHatDot_3, SEqHatDot_4; //estimated direction of the gyroscope error
    axis3f_t w_err;                                           //estimated direction of the gyroscope error (angular)
    axis3f_t h;                                               //computed flux in the earth frame

    // axulirary variables to avoid reapeated calculations
    float halfSEq_1 = 0.5f * q.q0;
    float halfSEq_2 = 0.5f * q.q1;
    float halfSEq_3 = 0.5f * q.q2;
    float halfSEq_4 = 0.5f * q.q3;
    float twoSEq_1 = 2.0f * q.q0;
    float twoSEq_2 = 2.0f * q.q1;
    float twoSEq_3 = 2.0f * q.q2;
    float twoSEq_4 = 2.0f * q.q3;
    float twob_x = 2.0f * b_x;
    float twob_z = 2.0f * b_z;
    float twob_xSEq_1 = 2.0f * b_x * q.q0;
    float twob_xSEq_2 = 2.0f * b_x * q.q1;
    float twob_xSEq_3 = 2.0f * b_x * q.q2;
    float twob_xSEq_4 = 2.0f * b_x * q.q3;
    float twob_zSEq_1 = 2.0f * b_z * q.q0;
    float twob_zSEq_2 = 2.0f * b_z * q.q1;
    float twob_zSEq_3 = 2.0f * b_z * q.q2;
    float twob_zSEq_4 = 2.0f * b_z * q.q3;
    float SEq_1SEq_2;
    float SEq_1SEq_3 = q.q0 * q.q2;
    float SEq_1SEq_4;
    float SEq_2SEq_3;
    float SEq_2SEq_4 = q.q1 * q.q3;
    float SEq_3SEq_4;
    float SEq_1SEq_1;
    float SEq_2SEq_2;
    float SEq_3SEq_3;
    float SEq_4SEq_4;
    float twom_x = 2.0f * mag.x;
    float twom_y = 2.0f * mag.y;
    float twom_z = 2.0f * mag.z;

    // normalize the accelerometer measurement
    inv_norm = -fastInvSqrt(accel.x * accel.x + accel.y * accel.y + accel.z * accel.z);
    if (isnan(inv_norm) || isinf(inv_norm)) {
        inv_norm = -1.f;
    }
    accel.x *= inv_norm;
    accel.y *= inv_norm;
    accel.z *= inv_norm;

    // normalize the magnetometer measurement
    inv_norm = fastInvSqrt(mag.x * mag.x + mag.y * mag.y + mag.z * mag.z);
    if (isnan(inv_norm) || isinf(inv_norm)) {
        inv_norm = 1.f;
    }
    mag.x *= inv_norm;
    mag.y *= inv_norm;
    mag.z *= inv_norm;

    // compute the objective function and Jacobian
    f_1 = twoSEq_2 * q.q3 - twoSEq_1 * q.q2 - accel.x;
    f_2 = twoSEq_1 * q.q1 + twoSEq_3 * q.q3 - accel.y;
    f_3 = 1.0f - twoSEq_2 * q.q1 - twoSEq_3 * q.q2 - accel.z;
    f_4 = twob_x * (0.5f - q.q2 * q.q2 - q.q3 * q.q3) + twob_z * (SEq_2SEq_4 - SEq_1SEq_3) - mag.x;
    f_5 = twob_x * (q.q1 * q.q2 - q.q0 * q.q3) + twob_z * (q.q0 * q.q1 + q.q2 * q.q3) - mag.y;
    f_6 = twob_x * (SEq_1SEq_3 + SEq_2SEq_4) + twob_z * (0.5f - q.q1 * q.q1 - q.q2 * q.q2) - mag.z;
    J_11or24 = twoSEq_3;
    J_12or23 = 2.0f * q.q3;
    J_13or22 = twoSEq_1;
    J_14or21 = twoSEq_2;
    J_32 = 2.0f * J_14or21;
    J_33 = 2.0f * J_11or24;
    J_41 = twob_zSEq_3;
    J_42 = twob_zSEq_4;
    J_43 = 2.0f * twob_xSEq_3 + twob_zSEq_1;
    J_44 = 2.0f * twob_xSEq_4 - twob_zSEq_2;
    J_51 = twob_xSEq_4 - twob_zSEq_2;
    J_52 = twob_xSEq_3 + twob_zSEq_1;
    J_53 = twob_xSEq_2 + twob_zSEq_4;
    J_54 = twob_xSEq_1 - twob_zSEq_3;
    J_61 = twob_xSEq_3;
    J_62 = twob_xSEq_4 - 2.0f * twob_zSEq_2;
    J_63 = twob_xSEq_1 - 2.0f * twob_zSEq_3;
    J_64 = twob_xSEq_2;

    // compute the gradient (matrix multiplication)
    SEqHatDot_1 = J_14or21 * f_2 - J_11or24 * f_1 - J_41 * f_4 - J_51 * f_5 + J_61 * f_6;
    SEqHatDot_2 = J_12or23 * f_1 + J_13or22 * f_2 - J_32 * f_3 + J_42 * f_4 + J_52 * f_5 + J_62 * f_6;
    SEqHatDot_3 = J_12or23 * f_2 - J_33 * f_3 - J_13or22 * f_1 - J_43 * f_4 + J_53 * f_5 + J_63 * f_6;
    SEqHatDot_4 = J_14or21 * f_1 + J_11or24 * f_2 - J_44 * f_4 - J_54 * f_5 + J_64 * f_6;

    // normalise the gradient to estimate direction of the gyroscope error
    inv_norm = fastInvSqrt(SEqHatDot_1 * SEqHatDot_1 + SEqHatDot_2 * SEqHatDot_2 + SEqHatDot_3 * SEqHatDot_3
                           + SEqHatDot_4 * SEqHatDot_4);
    if (isnan(inv_norm) || isinf(inv_norm)) {
        inv_norm = 1.f;
    }
    SEqHatDot_1 = SEqHatDot_1 * inv_norm;
    SEqHatDot_2 = SEqHatDot_2 * inv_norm;
    SEqHatDot_3 = SEqHatDot_3 * inv_norm;
    SEqHatDot_4 = SEqHatDot_4 * inv_norm;

    // compute angular estimated direction of the gyroscope error
    w_err.x = twoSEq_1 * SEqHatDot_2 - twoSEq_2 * SEqHatDot_1 - twoSEq_3 * SEqHatDot_4 + twoSEq_4 * SEqHatDot_3;
    w_err.y = twoSEq_1 * SEqHatDot_3 + twoSEq_2 * SEqHatDot_4 - twoSEq_3 * SEqHatDot_1 - twoSEq_4 * SEqHatDot_2;
    w_err.z = twoSEq_1 * SEqHatDot_4 - twoSEq_2 * SEqHatDot_3 + twoSEq_3 * SEqHatDot_2 - twoSEq_4 * SEqHatDot_1;

    // compute and remove the gyroscope biases
    w_b.x += w_err.x * configAHRS_MADGWICK_LOOP_TIME_S * configAHRS_MADGWICK_GYRO_DRIFT;
    w_b.y += w_err.y * configAHRS_MADGWICK_LOOP_TIME_S * configAHRS_MADGWICK_GYRO_DRIFT;
    w_b.z += w_err.z * configAHRS_MADGWICK_LOOP_TIME_S * configAHRS_MADGWICK_GYRO_DRIFT;
    gyro.x -= w_b.x;
    gyro.y -= w_b.y;
    gyro.z -= w_b.z;

    // compute the quaternion rate measured by gyroscopes
    SEqDot_omega_1 = -halfSEq_2 * gyro.x - halfSEq_3 * gyro.y - halfSEq_4 * gyro.z;
    SEqDot_omega_2 = halfSEq_1 * gyro.x + halfSEq_3 * gyro.z - halfSEq_4 * gyro.y;
    SEqDot_omega_3 = halfSEq_1 * gyro.y - halfSEq_2 * gyro.z + halfSEq_4 * gyro.x;
    SEqDot_omega_4 = halfSEq_1 * gyro.z + halfSEq_2 * gyro.y - halfSEq_3 * gyro.x;

    // compute then integrate the estimated quaternion rate
    q.q0 += (SEqDot_omega_1 - (configAHRS_MADGWICK_GYRO_ERROR * SEqHatDot_1)) * configAHRS_MADGWICK_LOOP_TIME_S;
    q.q1 += (SEqDot_omega_2 - (configAHRS_MADGWICK_GYRO_ERROR * SEqHatDot_2)) * configAHRS_MADGWICK_LOOP_TIME_S;
    q.q2 += (SEqDot_omega_3 - (configAHRS_MADGWICK_GYRO_ERROR * SEqHatDot_3)) * configAHRS_MADGWICK_LOOP_TIME_S;
    q.q3 += (SEqDot_omega_4 - (configAHRS_MADGWICK_GYRO_ERROR * SEqHatDot_4)) * configAHRS_MADGWICK_LOOP_TIME_S;

    /* Normalise quaternion */
    quaternionNorm(&q);

    /* Compute flux in the earth frame */
    SEq_1SEq_2 = q.q0 * q.q1;
    SEq_1SEq_3 = q.q0 * q.q2;
    SEq_1SEq_4 = q.q0 * q.q3;
    SEq_3SEq_4 = q.q2 * q.q3;
    SEq_2SEq_3 = q.q1 * q.q2;
    SEq_2SEq_4 = q.q1 * q.q3;
    SEq_1SEq_1 = q.q0 * q.q0;
    SEq_2SEq_2 = q.q1 * q.q1;
    SEq_3SEq_3 = q.q2 * q.q2;
    SEq_4SEq_4 = q.q3 * q.q3;
    h.x = twom_x * (0.5f - SEq_3SEq_3 - SEq_4SEq_4) + twom_y * (SEq_2SEq_3 - SEq_1SEq_4)
          + twom_z * (SEq_2SEq_4 + SEq_1SEq_3);
    h.y = twom_x * (SEq_2SEq_3 + SEq_1SEq_4) + twom_y * (0.5f - SEq_2SEq_2 - SEq_4SEq_4)
          + twom_z * (SEq_3SEq_4 - SEq_1SEq_2);
    h.z = twom_x * (SEq_2SEq_4 - SEq_1SEq_3) + twom_y * (SEq_3SEq_4 + SEq_1SEq_2)
          + twom_z * (0.5f - SEq_2SEq_2 - SEq_3SEq_3);

    /* Normalize the flux vector to have only components in the x and z */
    b_x = sqrtf((h.x * h.x) + (h.y * h.y));
    b_z = h.z;

    /* Convert quaterionion to Euler angles */
    angles->x = atan2f((2.f * (SEq_3SEq_4 + SEq_1SEq_2)), (SEq_1SEq_1 - SEq_2SEq_2 - SEq_3SEq_3 + SEq_4SEq_4));
    angles->y = -asinf((2.f * (SEq_2SEq_4 - SEq_1SEq_3)));
    angles->z = atan2f((2.f * (SEq_2SEq_3 + SEq_1SEq_4)), (SEq_1SEq_1 + SEq_2SEq_2 - SEq_3SEq_3 - SEq_4SEq_4));

    return;
}