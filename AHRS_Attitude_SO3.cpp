/****************************************************************************
 *
 *   Copyright (c) 2013 PX4 Development Team. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 * 3. Neither the name PX4 nor the names of its contributors ma_y be
 *    used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY Wa_y OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 ****************************************************************************/

/*
 * @file attitude_estimator_so3_main.cpp
 *
 * @author Hyon Lim <limhyon@gmail.com>
 * @author Anton Babushkin <anton.babushkin@me.com>
 *
 * Implementation of nonlinear complementary filters on the SO(3).
 * This code performs attitude estimation by using accelerometer, g_yroscopes and magnetometer.
 * Result is provided as quaternion, 1-2-3 Euler angle and rotation matrix.
 *
 * Theory of nonlinear complementary filters on the SO(3) is based on [1].
 * Quaternion realization of [1] is based on [2].
 * Optmized quaternion update code is based on Sebastian Madgwick's implementation.
 *
 * References
 *  [1] Mahony, R.; Hamel, T.; Pflimlin, Jean-Michel, "Nonlinear Complementary Filters on the Special Orthogonal Group," Automatic Control, IEEE Transactions on , vol.53, no.5, pp.1203,1218, June 2008
 *  [2] Euston, M.; Coote, P.; Mahony, R.; Jonghyuk Kim; Hamel, T., "A complementary filter for attitude estimation of a fixed-wing UAV," Intelligent Robots and Systems, 2008. IROS 2008. IEEE/RSJ International Conference on , vol., no., pp.340,345, 22-26 Sept. 2008
 */


#include "AHRS_Attitude_SO3.h"

//=====================================Constructor===========================================//
AHRS_Attitude_SO3::AHRS_Attitude_SO3(float loop_time_s, float Kp, float Ki){
    
    // Copy variables
    _Kp = Kp;
    _Ki = Ki;
    _loop_time_s = loop_time_s;
    
    // Initialize quaternion
    q0 = 1;
    q1 = 0;
    q2 = 0;
    q3 = 0;
    
    // Auxiliary variables to avoid repeated arithmetic
    _q0q0 = q0 * q0;
    _q0q1 = q0 * q1;
    _q0q2 = q0 * q2;
    _q0q3 = q0 * q3;
    _q1q1 = q1 * q1;
    _q1q2 = q1 * q2;
   	_q1q3 = q1 * q3;
    _q2q2 = q2 * q2;
    _q2q3 = q2 * q3;
    _q3q3 = q3 * q3;
    return;
}

//====================================Public Members==========================================//
//--------------------------Compute---------------------------//
void AHRS_Attitude_SO3::compute(float g_x, float g_y, float g_z, float a_x, float a_y, float a_z, float m_x, float m_y, float m_z){
    
    float inv_norm;
    float halfex = 0.0f, halfey = 0.0f, halfez = 0.0f;
    
    //! If magnetometer measurement is available, use it.
    if(!((m_x == 0.0f) && (m_y == 0.0f) && (m_z == 0.0f))) {
        float hx, hy, hz, bx, bz;
        float halfwx, halfwy, halfwz;
        
        // Normalise magnetometer measurement
        // Will sqrt work better? PX4 system is powerful enough?
        inv_norm = 1.f / sqrtf(m_x * m_x + m_y * m_y + m_z * m_z);
        m_x *= inv_norm;
        m_y *= inv_norm;
        m_z *= inv_norm;
        
        // Reference direction of Earth's magnetic field
        hx = 2.0f * (m_x * (0.5f - _q2q2 - _q3q3) + m_y * (_q1q2 - _q0q3) + m_z * (_q1q3 + _q0q2));
        hy = 2.0f * (m_x * (_q1q2 + _q0q3) + m_y * (0.5f - _q1q1 - _q3q3) + m_z * (_q2q3 - _q0q1));
        hz = 2.0f * m_x * (_q1q3 - _q0q2) + 2.0f * m_y * (_q2q3 + _q0q1) + 2.0f * m_z * (0.5f - _q1q1 - _q2q2);
        bx = sqrtf(hx * hx + hy * hy);
        bz = hz;
        
        // Estimated direction of magnetic field
        halfwx = bx * (0.5f - _q2q2 - _q3q3) + bz * (_q1q3 - _q0q2);
        halfwy = bx * (_q1q2 - _q0q3) + bz * (_q0q1 + _q2q3);
        halfwz = bx * (_q0q2 + _q1q3) + bz * (0.5f - _q1q1 - _q2q2);
        
        // Error is sum of cross product between estimated direction and measured direction of field vectors
        halfex += (m_y * halfwz - m_z * halfwy);
        halfey += (m_z * halfwx - m_x * halfwz);
        halfez += (m_x * halfwy - m_y * halfwx);
    }
    
    // Compute feedback only if accelerometer measurement valid (avoids NaN in accelerometer normalisation)
    if(!((a_x == 0.0f) && (a_y == 0.0f) && (a_z == 0.0f))) {
        float halfvx, halfvy, halfvz;
        
        // Normalise accelerometer measurement
        inv_norm = 1.f / sqrtf(a_x * a_x + a_y * a_y + a_z * a_z);
        a_x *= inv_norm;
        a_y *= inv_norm;
        a_z *= inv_norm;
        
        // Estimated direction of gravity and magnetic field
        halfvx = _q1q3 - _q0q2;
        halfvy = _q0q1 + _q2q3;
        halfvz = _q0q0 - 0.5f + _q3q3;
        
        // Error is sum of cross product between estimated direction and measured direction of field vectors
        halfex += a_y * halfvz - a_z * halfvy;
        halfey += a_z * halfvx - a_x * halfvz;
        halfez += a_x * halfvy - a_y * halfvx;
    }
    
    // Apply feedback only when valid data has been gathered from the accelerometer or magnetometer
    if(halfex != 0.0f && halfey != 0.0f && halfez != 0.0f) {
        // Compute and apply integral feedback if enabled
        if(_Ki > 0.0f) {
            _gyro_bias[0] += _Ki * halfex * _loop_time_s;	// integral error scaled by Ki
            _gyro_bias[1] += _Ki * halfey * _loop_time_s;
            _gyro_bias[2] += _Ki * halfez * _loop_time_s;
            
            // apply integral feedback
            g_x += _gyro_bias[0];
            g_y += _gyro_bias[1];
            g_z += _gyro_bias[2];
        }
        else {
            _gyro_bias[0] = 0.0f;	// prevent integral windup
            _gyro_bias[1] = 0.0f;
            _gyro_bias[2] = 0.0f;
        }
        
        // Apply proportional feedback
        g_x += _Kp * halfex;
        g_y += _Kp * halfey;
        g_z += _Kp * halfez;
    }
    
    //! Integrate rate of change of quaternion
#if 0
    g_x *= (0.5f * _loop_time_s);		// pre-multiply common factors
    g_y *= (0.5f * _loop_time_s);
    g_z *= (0.5f * _loop_time_s);
#endif
    
    // Time derivative of quaternion. q_dot = 0.5*q\otimes omega.
    //! q_k = q_{k-1} + _loop_time_s*\dot{q}
    //! \dot{q} = 0.5*q \otimes P(\omega)
    float dq0 = 0.5f * (- q1 * g_x - q2 * g_y - q3 * g_z);
    float dq1 = 0.5f * (q0 * g_x + q2 * g_z - q3 * g_y);
    float dq2 = 0.5f * (q0 * g_y - q1 * g_z + q3 * g_x);
    float dq3 = 0.5f * (q0 * g_z + q1 * g_y - q2 * g_x);
    
    q0 += _loop_time_s * dq0;
    q1 += _loop_time_s * dq1;
    q2 += _loop_time_s * dq2;
    q3 += _loop_time_s * dq3;
    
    // Normalise quaternion
    inv_norm = 1.f / sqrtf(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
    q0 *= inv_norm;
    q1 *= inv_norm;
    q2 *= inv_norm;
    q3 *= inv_norm;
    
    // Auxiliary variables to avoid repeated arithmetic
    _q0q0 = q0 * q0;
    _q0q1 = q0 * q1;
    _q0q2 = q0 * q2;
    _q0q3 = q0 * q3;
    _q1q1 = q1 * q1;
    _q1q2 = q1 * q2;
   	_q1q3 = q1 * q3;
    _q2q2 = q2 * q2;
    _q2q3 = q2 * q3;
    _q3q3 = q3 * q3;
    
    // Convert from quaternion to Euler angles
    Roll = atan2f((2.f * (_q2q3 + _q0q1)), (_q0q0 - _q1q1 - _q2q2 + _q3q3));
    Pitch = -asinf((2.f * (_q1q3 - _q0q2)));
    Yaw = atan2f((2.f * (_q1q2 + _q0q3)), (_q0q0 + _q1q1 - _q2q2 - _q3q3));
    return;

}

//-----------------------------------Starting values-----------------------------------------//
void AHRS_Attitude_SO3::set_starting_values(float Phi0, float Theta0, float Psi0){
    float cphi = cosf(Phi0 * 0.5f);
    float sphi = sinf(Phi0 * 0.5f);
    
    float ctheta = cosf(Theta0 * 0.5f);
    float stheta = sinf(Theta0 * 0.5f);
    
    float cpsi = cosf(Psi0 * 0.5f);
    float spsi = sinf(Psi0 * 0.5f);
    
    q0 = cphi * ctheta * cpsi + sphi * stheta * spsi;
    q1 = sphi * ctheta * cpsi - cphi * stheta * spsi;
    q2 = cphi * stheta * cpsi + sphi * ctheta * spsi;
    q3 = cphi * ctheta * spsi - sphi * stheta * cpsi;
    
    // Auxiliary variables to avoid repeated arithmetic
    _q0q0 = q0 * q0;
    _q0q1 = q0 * q1;
    _q0q2 = q0 * q2;
    _q0q3 = q0 * q3;
    _q1q1 = q1 * q1;
    _q1q2 = q1 * q2;
    _q1q3 = q1 * q3;
    _q2q2 = q2 * q2;
    _q2q3 = q2 * q3;
    _q3q3 = q3 * q3;
    return;

}