//  IMU_Madgwick.cpp
//
//
//  Created by Sebastian Madgwick
//
//
#include "IMU_Madgwick.h"

//=====================================Constructor===========================================//
IMU_Madgwick::IMU_Madgwick(float gyro_error) {
    // Initial quaternion
    q1 = 1f;
    q2 = 0f;
    q3 = 0f;
    q4 = 0f;
    // filter parameters
    _beta = sqrtf(3.0f * 0.25f) * gyro_error;
}

//====================================Public Members==========================================//
//--------------------------Compute---------------------------//
void IMU_Madgwick::compute(float g_x, float g_y, float g_z, float a_x, float a_y, float a_z, float loop_time_s) {
    // Local system variables
    float inv_norm; //vector norm
    float SEqDot_omega_1, SEqDot_omega_2, SEqDot_omega_3,
        SEqDot_omega_4;                                       //quaternion derivative from gyroscope elements
    float f_1, f_2, f_3;                                      //objective function elements
    float J_11or24, J_12or23, J_13or22, J_14or21, J_32, J_33; //objective function Jacobian elements
    float SEqHatDot_1, SEqHatDot_2, SEqHatDot_3, SEqHatDot_4; //estimated direction of the gyroscope error

    // axulirary variables to avoid reapeated calculations
    float halfSEq_1 = 0.5f * q1;
    float halfSEq_2 = 0.5f * q2;
    float halfSEq_3 = 0.5f * q3;
    float halfSEq_4 = 0.5f * q4;
    float twoSEq_1 = 2.0f * q1;
    float twoSEq_2 = 2.0f * q2;
    float twoSEq_3 = 2.0f * q3;

    /* Re-orientate readings */
    float tmp = a_x;
    a_x = a_y;
    a_y = tmp;
    a_z *= -1;
    tmp = g_x;
    g_x = g_y;
    g_y = tmp;
    g_z *= -1;

    // normalize the accelerometer measurement
    inv_norm = 1.f / sqrtf(a_x * a_x + a_y * a_y + a_z * a_z);
    if (isnan(inv_norm) || isinf(inv_norm)) {
        inv_norm = 1.f;
    }
    a_x *= inv_norm;
    a_y *= inv_norm;
    a_z *= inv_norm;

    // compute the objective function and Jacobian
    f_1 = twoSEq_2 * q4 - twoSEq_1 * q3 - a_x;
    f_2 = twoSEq_1 * q2 + twoSEq_3 * q4 - a_y;
    f_3 = 1.0f - twoSEq_2 * q2 - twoSEq_3 * q3 - a_z;
    J_11or24 = twoSEq_3;
    J_12or23 = 2.0f * q4;
    J_13or22 = twoSEq_1;
    J_14or21 = twoSEq_2;
    J_32 = 2.0f * J_14or21;
    J_33 = 2.0f * J_11or24;

    // compute the gradient (matrix multiplication)
    SEqHatDot_1 = J_14or21 * f_2 - J_11or24 * f_1;
    SEqHatDot_2 = J_12or23 * f_1 + J_13or22 * f_2 - J_32 * f_3;
    SEqHatDot_3 = J_12or23 * f_2 - J_33 * f_3 - J_13or22 * f_1;
    SEqHatDot_4 = J_14or21 * f_1 + J_11or24 * f_2;

    // normalise the gradient to estimate direction of the gyroscope error
    inv_norm = 1.f
               / sqrtf(SEqHatDot_1 * SEqHatDot_1 + SEqHatDot_2 * SEqHatDot_2 + SEqHatDot_3 * SEqHatDot_3
                       + SEqHatDot_4 * SEqHatDot_4);
    if (isnan(inv_norm) || isinf(inv_norm)) {
        inv_norm = 1.f;
    }
    SEqHatDot_1 *= inv_norm;
    SEqHatDot_2 *= inv_norm;
    SEqHatDot_3 *= inv_norm;
    SEqHatDot_4 *= inv_norm;

    // compute the quaternion rate measured by gyroscopes
    SEqDot_omega_1 = -halfSEq_2 * g_x - halfSEq_3 * g_y - halfSEq_4 * g_z;
    SEqDot_omega_2 = halfSEq_1 * g_x + halfSEq_3 * g_z - halfSEq_4 * g_y;
    SEqDot_omega_3 = halfSEq_1 * g_y - halfSEq_2 * g_z + halfSEq_4 * g_x;
    SEqDot_omega_4 = halfSEq_1 * g_z + halfSEq_2 * g_y - halfSEq_3 * g_x;

    // compute then integrate the estimated quaternion rate
    q1 += (SEqDot_omega_1 - (_beta * SEqHatDot_1)) * loop_time_s;
    q2 += (SEqDot_omega_2 - (_beta * SEqHatDot_2)) * loop_time_s;
    q3 += (SEqDot_omega_3 - (_beta * SEqHatDot_3)) * loop_time_s;
    q4 += (SEqDot_omega_4 - (_beta * SEqHatDot_4)) * loop_time_s;

    // normalise quaternion
    inv_norm = 1.f / sqrtf(q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4);
    if (isnan(inv_norm) || isinf(inv_norm)) {
        inv_norm = 1.f;
    }
    q1 *= inv_norm;
    q2 *= inv_norm;
    q3 *= inv_norm;
    q4 *= inv_norm;

    // convert from quaternion to Euler angles
    // switched x and y and changed sign to z to align with current reference system
    Pitch = atan2f((2.f * (q3 * q4 + q1 * q2)), (q1 * q1 - q2 * q2 - q3 * q3 + q4 * q4));
    Roll = -asinf((2.f * (q2 * q4 - q1 * q3)));
    return;
}