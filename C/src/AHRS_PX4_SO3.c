/* BEGIN Header */
/**
 ******************************************************************************
 * \file            AHRS_PX4_SO3.c
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

/* Includes ------------------------------------------------------------------*/

#include "AHRS_PX4_SO3.h"
#include "math.h"
#include "quaternion.h"

/* Private variables ---------------------------------------------------------*/
static quaternion_t _q = {1, 0, 0, 0};
static float _q0q0 = 1.f;
static float _q0q1 = 0.f;
static float _q0q2 = 0.f;
static float _q0q3 = 0.f;
static float _q1q1 = 0.f;
static float _q1q2 = 0.f;
static float _q1q3 = 0.f;
static float _q2q2 = 0.f;
static float _q2q3 = 0.f;
static float _q3q3 = 0.f;

/* Functions -----------------------------------------------------------------*/
void AHRS_PX4_S03_update(axis3f_t* angles, axis3f_t accel, axis3f_t gyro, axis3f_t mag) {
    static axis3f_t gyro_bias = {0, 0, 0};
    float inv_norm;
    axis3f_t halfe = {0, 0, 0};

    /* If magnetometer measurement is available, use it */
    if (!((mag.x == 0.0f) && (mag.y == 0.0f) && (mag.z == 0.0f))) {
        axis3f_t h, halfw;
        float bx, bz;

        /* Normalise magnetometer measurement */
        /* Will sqrt work better? PX4 system is powerful enough? */
        inv_norm = 1.f / sqrtf(mag.x * mag.x + mag.y * mag.y + mag.z * mag.z);
        mag.x *= inv_norm;
        mag.y *= inv_norm;
        mag.z *= inv_norm;

        /* Reference direction of Earth's magnetic field */
        h.x = 2.0f * (mag.x * (0.5f - _q2q2 - _q3q3) + mag.y * (_q1q2 - _q0q3) + mag.z * (_q1q3 + _q0q2));
        h.y = 2.0f * (mag.x * (_q1q2 + _q0q3) + mag.y * (0.5f - _q1q1 - _q3q3) + mag.z * (_q2q3 - _q0q1));
        h.z = 2.0f * mag.x * (_q1q3 - _q0q2) + 2.0f * mag.y * (_q2q3 + _q0q1) + 2.0f * mag.z * (0.5f - _q1q1 - _q2q2);
        bx = sqrtf(h.x * h.x + h.y * h.y);
        bz = h.z;

        /* Estimated direction of magnetic field */
        halfw.x = bx * (0.5f - _q2q2 - _q3q3) + bz * (_q1q3 - _q0q2);
        halfw.y = bx * (_q1q2 - _q0q3) + bz * (_q0q1 + _q2q3);
        halfw.z = bx * (_q0q2 + _q1q3) + bz * (0.5f - _q1q1 - _q2q2);

        /* Error is sum of cross product between estimated direction and measured direction of field vectors */
        halfe.x += (mag.y * halfw.z - mag.z * halfw.y);
        halfe.y += (mag.z * halfw.x - mag.x * halfw.z);
        halfe.z += (mag.x * halfw.y - mag.y * halfw.x);
    }

    /* Compute feedback only if accelerometer measurement valid (avoids NaN in accelerometer normalisation) */
    if (!((accel.x == 0.0f) && (accel.y == 0.0f) && (accel.z == 0.0f))) {
        axis3f_t halfv;

        /* Normalise accelerometer measurement */
        inv_norm = 1.f / sqrtf(accel.x * accel.x + accel.y * accel.y + accel.z * accel.z);
        accel.x *= inv_norm;
        accel.y *= inv_norm;
        accel.z *= inv_norm;

        /* Estimated direction of gravity and magnetic field */
        halfv.x = _q1q3 - _q0q2;
        halfv.y = _q0q1 + _q2q3;
        halfv.z = _q0q0 - 0.5f + _q3q3;

        /* Error is sum of cross product between estimated direction and measured direction of field vectors */
        halfe.x += accel.y * halfv.z - accel.z * halfv.y;
        halfe.y += accel.z * halfv.x - accel.x * halfv.z;
        halfe.z += accel.x * halfv.y - accel.y * halfv.x;
    }

    /* Apply feedback only when valid data has been gathered from the accelerometer or magnetometer */
    if (halfe.x != 0.0f && halfe.y != 0.0f && halfe.z != 0.0f) {
        /* Compute and apply integral feedback if enabled */
        if (configAHRS_PX4_SO3_KI > 0) {
            /* Integral error scaled by Ki */
            gyro_bias.x += configAHRS_PX4_SO3_KI * halfe.x * configAHRS_PX4_SO3_LOOP_TIME_S;
            gyro_bias.y += configAHRS_PX4_SO3_KI * halfe.y * configAHRS_PX4_SO3_LOOP_TIME_S;
            gyro_bias.z += configAHRS_PX4_SO3_KI * halfe.z * configAHRS_PX4_SO3_LOOP_TIME_S;

            /* Apply integral feedback */
            gyro.x += gyro_bias.x;
            gyro.y += gyro_bias.y;
            gyro.z += gyro_bias.z;
        } else {
            gyro_bias.x = 0.0f; // prevent integral windup
            gyro_bias.y = 0.0f;
            gyro_bias.z = 0.0f;
        }

        /* Apply proportional feedback */
        gyro.x += configAHRS_PX4_SO3_KP * halfe.x;
        gyro.y += configAHRS_PX4_SO3_KP * halfe.y;
        gyro.z += configAHRS_PX4_SO3_KP * halfe.z;
    }

    /* Integrate rate of change of quaternion */
#if 0
    gyro.x *= (0.5f * configAHRS_PX4_SO3_LOOP_TIME_S);		// pre-multiply common factors
    gyro.y *= (0.5f * configAHRS_PX4_SO3_LOOP_TIME_S);
    gyro.z *= (0.5f * configAHRS_PX4_SO3_LOOP_TIME_S);
#endif

    /* Time derivative of quaternion. q_dot = 0.5*q\otimes omega */
    /* q_k = q_{k-1} + configAHRS_PX4_SO3_LOOP_TIME_S*\dot{q} */
    /* \dot{q} = 0.5*q \otimes P(\omega) */
    float dq0 = 0.5f * (-_q.q1 * gyro.x - _q.q2 * gyro.y - _q.q3 * gyro.z);
    float dq1 = 0.5f * (_q.q0 * gyro.x + _q.q2 * gyro.z - _q.q3 * gyro.y);
    float dq2 = 0.5f * (_q.q0 * gyro.y - _q.q1 * gyro.z + _q.q3 * gyro.x);
    float dq3 = 0.5f * (_q.q0 * gyro.z + _q.q1 * gyro.y - _q.q2 * gyro.x);

    _q.q0 += configAHRS_PX4_SO3_LOOP_TIME_S * dq0;
    _q.q1 += configAHRS_PX4_SO3_LOOP_TIME_S * dq1;
    _q.q2 += configAHRS_PX4_SO3_LOOP_TIME_S * dq2;
    _q.q3 += configAHRS_PX4_SO3_LOOP_TIME_S * dq3;

    /* Normalise quaternion */
    quaternionNorm(&_q);

    /* Auxiliary variables to avoid repeated arithmetic */
    _q0q0 = _q.q0 * _q.q0;
    _q0q1 = _q.q0 * _q.q1;
    _q0q2 = _q.q0 * _q.q2;
    _q0q3 = _q.q0 * _q.q3;
    _q1q1 = _q.q1 * _q.q1;
    _q1q2 = _q.q1 * _q.q2;
    _q1q3 = _q.q1 * _q.q3;
    _q2q2 = _q.q2 * _q.q2;
    _q2q3 = _q.q2 * _q.q3;
    _q3q3 = _q.q3 * _q.q3;

    /* Convert quaterionion to Euler angles */
    quaternionToEuler(&_q, angles);

    return;
}

void AHRS_PX4_S03_reset(axis3f_t* angles, float phi0, float theta0, float psi0) {
    float cPhi = cosf(phi0 * 0.5f);
    float sPhi = sinf(phi0 * 0.5f);

    float cTheta = cosf(theta0 * 0.5f);
    float sTheta = sinf(theta0 * 0.5f);

    float cPsi = cosf(psi0 * 0.5f);
    float sPsi = sinf(psi0 * 0.5f);

    _q.q0 = cPhi * cTheta * cPsi + sPhi * sTheta * sPsi;
    _q.q1 = sPhi * cTheta * cPsi - cPhi * sTheta * sPsi;
    _q.q2 = cPhi * sTheta * cPsi + sPhi * cTheta * sPsi;
    _q.q3 = cPhi * cTheta * sPsi - sPhi * sTheta * cPsi;

    /* Auxiliary variables to avoid repeated arithmetic */
    _q0q0 = _q.q0 * _q.q0;
    _q0q1 = _q.q0 * _q.q1;
    _q0q2 = _q.q0 * _q.q2;
    _q0q3 = _q.q0 * _q.q3;
    _q1q1 = _q.q1 * _q.q1;
    _q1q2 = _q.q1 * _q.q2;
    _q1q3 = _q.q1 * _q.q3;
    _q2q2 = _q.q2 * _q.q2;
    _q2q3 = _q.q2 * _q.q3;
    _q3q3 = _q.q3 * _q.q3;

    /* Convert quaterionion to Euler angles */
    quaternionToEuler(&_q, angles);

    return;
}