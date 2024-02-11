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
 * 3. Neither the name PX4 nor the names of its contributors may be
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
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
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
 * This code performs attitude estimation by using accelerometer, gyroscopes and magnetometer.
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


#ifndef _AHRS_ATTITUDE_SO3_H_
#define _AHRS_ATTITUDE_SO3_H_
#include "math.h"

//======================================Parameters=============================================//
#define AHRS_Attitude_SO3_Kp    1.f     //default Kp value
#define AHRS_Attitude_SO3_Ki    0.05f   //default Ki value

class AHRS_Attitude_SO3 {
public:
    float Roll, Pitch, Yaw; //estimated orientation (Euler angles in RPY order: Phi=roll, Teta=pitch, Psi=yaw), in rad
    float q0, q1, q2, q3; //estimated quaternion
    AHRS_Attitude_SO3(float loop_time_s, float Kp = AHRS_Attitude_SO3_Kp, float Ki = AHRS_Attitude_SO3_Ki); //constructor
    void compute(float g_x, float g_y, float g_z, float a_x, float a_y, float a_z, float m_x, float m_y, float m_z);
    void set_starting_values(float Phi0, float Theta0, float Psi0);
    
private:
    float _gyro_bias[3]={0.f, 0.f, 0.f}; //gyro bias vector
    float _Kp, _Ki; //gain values
    float _loop_time_s; //loop time, in s
    float _q0q0, _q0q1, _q0q2, _q0q3, _q1q1, _q1q2, _q1q3, _q2q2, _q2q3, _q3q3; //auxiliary variables
};

#endif
