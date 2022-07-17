/****************************************************************************
 *
 *   Copyright (c) 2012-2014 PX4 Development Team. All rights reserved.
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
 * ANY u(4,0) OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 ****************************************************************************/

/*
 * @file AHRS_Attitude_EKF.cpp
 *
 * Extended Kalman Filter for Attitude Estimation.
 *
 * @author Tobias Naegeli <naegelit@student.ethz.ch>
 * @author Lorenz Meier <lm@inf.ethz.ch>
 * @author Thomas Gubler <thomasgubler@gmail.com>
 */


#include "AHRS_Attitude_EKF.h"

//=====================================Constructors==========================================//
AHRS_Attitude_EKF::AHRS_Attitude_EKF(float g_val, float loop_time_s):
u(12,1), _A(12,12), _C(3,12), _J(3,3), _Ji(3,3), _Q(12,12), _M(1,1), _K(1,1), _R_gyro(3,3), _R_accel(3,3), _R_mag(3,3), _P(12,12){
    MatrixXf inertia_m(3,3);
    inertia_m.identity();
    initialize(g_val, loop_time_s, loop_time_s, loop_time_s, loop_time_s, inertia_m);
    return;
};

AHRS_Attitude_EKF::AHRS_Attitude_EKF(float g_val, float loop_time_s, MatrixXf &inertia_m):
u(12,1), _A(12,12), _C(3,12), _J(3,3), _Ji(3,3), _Q(12,12), _M(1,1), _K(1,1), _R_gyro(3,3), _R_accel(3,3), _R_mag(3,3), _P(12,12){
    initialize(g_val, loop_time_s, loop_time_s, loop_time_s, loop_time_s, inertia_m);
    return;
};

AHRS_Attitude_EKF::AHRS_Attitude_EKF(float g_val, float loop_time_pred_s, float loop_time_update_gyro_s, float loop_time_update_acc_s, float loop_time_update_mag_s, MatrixXf &inertia_m):
u(12,1), _A(12,12), _C(3,12), _J(3,3), _Ji(3,3), _Q(12,12), _M(1,1), _K(1,1), _R_gyro(3,3), _R_accel(3,3), _R_mag(3,3), _P(12,12){
    initialize(g_val, loop_time_pred_s, loop_time_update_gyro_s, loop_time_update_acc_s, loop_time_update_mag_s, inertia_m);
    return;
};

//====================================Public Members==========================================//
//-------------------------------------Prediction---------------------------------------------//
void AHRS_Attitude_EKF::prediction(){
    /*
     wx=  u(0);   % x  body angular rate
     wy=  u(1);   % y  body angular rate
     wz=  u(2);   % z  body angular rate
     
     wax=  u(3);  % x  body angular acceleration
     way=  u(4);  % y  body angular acceleration
     waz=  u(5);  % z  body angular acceleration
     
     zex=  u(6);  % x  component gravity vector
     zey=  u(7);  % y  component gravity vector
     zez=  u(8);  % z  component gravity vector
     
     mux=  u(9); % x  component magnetic field vector
     muy=  u(10); % y  component magnetic field vector
     muz=  u(11); % z  component magnetic field vector*/
    
    float w1 = _loop_time_pred_s * u(0,0);
    float w2 = _loop_time_pred_s * u(1,0);
    float w3 = _loop_time_pred_s * u(2,0);
    //A matrix
    _A(0,0) = 1.f;
    _A(0,3) = _loop_time_pred_s;
    _A(1,1) = 1.f;
    _A(1,4) = _loop_time_pred_s;
    _A(2,2) = 1.f;
    _A(2,5) = _loop_time_pred_s;
    _A(3,3) = 1.f;
    _A(4,4) = 1.f;
    _A(5,5) = 1.f;
    _A(6,1) = - _loop_time_pred_s * u(8,0);
    _A(6,2) = _loop_time_pred_s * u(7,0);
    _A(6,6) = 1.f;
    _A(6,7) = w3;
    _A(6,8) = - w2;
    _A(7,0) = _loop_time_pred_s * u(8,0);
    _A(7,2) = - _loop_time_pred_s * u(6,0);
    _A(7,6) = - w3;
    _A(7,7) = 1.f;
    _A(7,8) = w1;
    _A(8,0) = - _loop_time_pred_s * u(7,0);
    _A(8,1) = _loop_time_pred_s * u(6,0);
    _A(8,6) = w2;
    _A(8,7) = - w1;
    _A(8,8) = 1.f;
    _A(9,1) = - _loop_time_pred_s * u(11,0);
    _A(9,2) = _loop_time_pred_s * u(10,0);
    _A(9,9) = 1.f;
    _A(9,10) = w3;
    _A(9,11) = - w2;
    _A(10,0) = _loop_time_pred_s * u(11,0);
    _A(10,2) = - _loop_time_pred_s * u(9,0);
    _A(10,9) = - w3;
    _A(10,10) = 1.f;
    _A(10,11) = w1;
    _A(11,0) = - _loop_time_pred_s * u(10,0);
    _A(11,1) = _loop_time_pred_s * u(9,0);
    _A(11,9) = w2;
    _A(11,10) = - w1;
    _A(11,11) = 1.f;
    
    //Predicted P matrix
    
    _P = QuadProd(_A,_P) + _Q;
    
    
    //Compute the apriori state estimate from the previous aposteriori estimate
    
    //Body angular rates prediction
#ifdef AHRS_ATTITUDE_EKF_USE_COMPLETE_INERTIA_MATRIX
    
    float t1 = _J(0,1) * u(3,0) + _J(1,1) * u(4,0) + _J(1,2) * u(5,0);
    float t2 = _J(0,0) * u(3,0) + _J(0,1) * u(4,0) + _J(0,2) * u(5,0);
    float t3 = _J(0,2) * u(3,0) + _J(1,2) * u(4,0) + _J(2,2) * u(5,0);
    float tmp1 = u(3,0) * t1 - u(4,0) * t2;
    float tmp2 = u(3,0) * t3 - u(5,0) * t2;
    float tmp3 = u(4,0) * t3 - u(5,0) * t1;

    float wax = u(3,0) - _loop_time_pred_s * (_Ji(0,2) * tmp1 - _Ji(0,1) * tmp2 + _Ji(0,0) * tmp3);
    float way = u(4,0) - _loop_time_pred_s * (_Ji(1,2) * tmp1 - _Ji(1,1) * tmp2 + _Ji(0,1) * tmp3);
    float waz = u(5,0) - _loop_time_pred_s * (_Ji(2,2) * tmp1 - _Ji(1,2) * tmp2 + _Ji(0,2) * tmp3);
    
#elif AHRS_ATTITUDE_EKF_USE_DIAGONAL_INERTIA_MATRIX
    
    float wax = u(3,0) + (_loop_time_pred_s * u(4,0) * u(5,0) * (_J(1,1) - _J(2,2))) / _J(0,0);
    float way = u(4,0) - (_loop_time_pred_s * u(3,0) * u(5,0) * (_J(0,0) - _J(2,2))) / _J(1,1);
    float waz = u(5,0) + (_loop_time_pred_s * u(3,0) * u(4,0) * (_J(0,0) - _J(1,1))) / _J(2,2);
    
#else
    
    float wax = u(3,0);
    float way = u(4,0);
    float waz = u(5,0);
    
#endif

    float delta_u6 = _loop_time_pred_s * (u(2,0) * u(7,0) - u(1,0) * u(8,0));
    float delta_u7 = _loop_time_pred_s * (u(0,0) * u(8,0) - u(2,0) * u(6,0));
    float delta_u8 = _loop_time_pred_s * (u(1,0) * u(6,0) - u(0,0) * u(7,0));
    
    float delta_u9 = _loop_time_pred_s * (u(2,0) * u(10,0) - u(1,0) * u(11,0));
    float delta_u10 = _loop_time_pred_s * (u(0,0) * u(11,0) - u(2,0) * u(9,0));
    float delta_u11 = _loop_time_pred_s * (u(1,0) * u(9,0) - u(0,0) * u(10,0));
    
    //Predict state
    u(0,0) += _loop_time_pred_s * wax;
    u(1,0) += _loop_time_pred_s * way;
    u(2,0) += _loop_time_pred_s * waz;
    u(3,0) = wax;
    u(4,0) = way;
    u(5,0) = waz;
    u(6,0) += delta_u6;
    u(7,0) += delta_u7;
    u(8,0) += delta_u8;
    u(9,0) += delta_u9;
    u(10,0) += delta_u10;
    u(11,0) += delta_u11;
    return;
}

//---------------------------------Update with gyro-----------------------------------------//
void AHRS_Attitude_EKF::update_gyro(float g_x, float g_y, float g_z){
    
    //C matrix
    _C.zeros();
    _C(0,0) = 1.f;
    _C(1,1) = 1.f;
    _C(2,2) = 1.f;
    
    //Delta measures
    MatrixXf Delta_m(3,1);
    Delta_m(0,0) = g_x - u(0,0);
    Delta_m(1,0) = g_y - u(1,0);
    Delta_m(2,0) = g_z - u(2,0);

    //Gain matrix K
    _M = QuadProd(_C,_P) + _R_gyro;
    _K = _P * (~_C) * (!_M);
    //Correct state vector
    u += _K * Delta_m;

    //Updated P matrix
    _P -= _K * _C * _P;
    //_P=(_P+(~_P))*0.5; //guarantees P to be symmetric
    return;
}

//---------------------------------Update with accel----------------------------------------//
void AHRS_Attitude_EKF::update_accel(float a_x, float a_y, float a_z){
    
    //C matrix
    _C.zeros();
    _C(0,6) = 1.f;
    _C(1,7) = 1.f;
    _C(2,8) = 1.f;
    
    //Delta measures
    MatrixXf Delta_m(3,1);
    Delta_m(0,0) = a_x - u(6,0);
    Delta_m(1,0) = a_y - u(7,0);
    Delta_m(2,0) = a_z - u(8,0);
    //mprint(Delta_m);
    //Gain matrix K
    _M = QuadProd(_C,_P) + _R_accel;
    _K = _P * (~_C) * (!_M);
    //Correct state vector
    u += _K * Delta_m;
    
    //Updated P matrix
    _P -= _K * _C * _P;
    //_P=(_P+(~_P))*0.5; //guarantees P to be symmetric
    return;
}

//----------------------------------Update with mag-----------------------------------------//
void AHRS_Attitude_EKF::update_mag(float m_x, float m_y, float m_z){
    
    //C matrix
    _C.zeros();
    _C(0,9) = 1.f;
    _C(1,10) = 1.f;
    _C(2,11) = 1.f;
    
    //Delta measures
    MatrixXf Delta_m(3,1);
    Delta_m(0,0) = m_x - u(9,0);
    Delta_m(1,0) = m_y - u(10,0);
    Delta_m(2,0) = m_z - u(11,0);
    
    //Gain matrix K
    _M = QuadProd(_C,_P) + _R_mag;
    _K = _P * (~_C) * (!_M);
    //Correct state vector
    u += _K * Delta_m;
    
    //Updated P matrix
    _P -= _K * _C * _P;
    //_P=(_P+(~_P))*0.5; //guarantees P to be symmetric
    return;
}

//--------------------------------Compute euler angles--------------------------------------//
void AHRS_Attitude_EKF::euler_angles(){
    
    //Normalize accelerometer estimate
    float inv_norm = 1.f / sqrt(u(6,0) * u(6,0) + u(7,0) * u(7,0) + u(8,0) * u(8,0));
    if (isnan(inv_norm)||isinf(inv_norm)){
        inv_norm=1.f/_g_val;
    }
    float ax = u(6,0) * inv_norm;
    float ay = u(7,0) * inv_norm;
    float az = u(8,0) * inv_norm;
    
    //Normalize magnetometer estimate
    inv_norm = 1.f / sqrt(u(9,0) * u(9,0) + u(10,0) * u(10,0) + u(11,0) * u(11,0));
    if (isnan(inv_norm)||isinf(inv_norm)){
        inv_norm=1.f;
    }
    float mx = u(9,0) * inv_norm;
    float my = u(10,0) * inv_norm;
    float mz = u(11,0) * inv_norm;
    
    //Compute Euler angles
    Roll = atan2f(-ay,-az);
    Pitch = asinf(ax);
    float sphi = sin(Roll);
    float cphi = cosf(Roll);
    float steta = ax;
    float cteta = cosf(Pitch);
    float Yh = my * cphi  - mz * sphi ;
    float Xh = mx * cteta + (my * sphi + mz * cphi) * steta;
    Yaw = atan2(-Yh, Xh);
    return;
}

//--------------------------------Set inclination angle--------------------------------------//
void AHRS_Attitude_EKF::set_inclination(float incl_angle){
    
    u(9,0)=cosf(incl_angle);
    u(10,0)=0;
    u(11,0)=sinf(incl_angle);
}

//----------------------------------Set process noises---------------------------------------//
void AHRS_Attitude_EKF::set_process_noise(float rot_sp, float rot_acc, float acc, float mag){
    
    _Q(0,0) = rot_sp * _loop_time_pred_s;
    _Q(1,1) = rot_sp * _loop_time_pred_s;
    _Q(2,2) = rot_sp * _loop_time_pred_s;
    _Q(3,3) = rot_acc * _loop_time_pred_s;
    _Q(4,4) = rot_acc * _loop_time_pred_s;
    _Q(5,5) = rot_acc * _loop_time_pred_s;
    _Q(6,6) = acc * _loop_time_pred_s;
    _Q(7,7) = acc * _loop_time_pred_s;
    _Q(8,8) = acc * _loop_time_pred_s;
    _Q(9,9) = mag * _loop_time_pred_s;
    _Q(10,10) = mag * _loop_time_pred_s;
    _Q(11,11) = mag * _loop_time_pred_s;
    return;
}

//------------------------------------Set gyro noises----------------------------------------//
void AHRS_Attitude_EKF::set_gyro_noise(float g){
    
    float inv_loop_time = 1.f / _loop_time_update_gyro_s;
    _R_gyro(0,0) = g * inv_loop_time;
    _R_gyro(1,1) = g * inv_loop_time;
    _R_gyro(2,2) = g * inv_loop_time;
    return;
}

//-------------------------------------Set acc noises----------------------------------------//
void AHRS_Attitude_EKF::set_accel_noise(float a){
    
    float inv_loop_time = 1.f / _loop_time_update_acc_s;
    _R_accel(0,0) = a * inv_loop_time;
    _R_accel(1,1) = a * inv_loop_time;
    _R_accel(2,2) = a * inv_loop_time;
    return;
}

//-------------------------------------Set mag noises----------------------------------------//
void AHRS_Attitude_EKF::set_mag_noise(float m){
    
    float inv_loop_time = 1.f / _loop_time_update_mag_s;
    _R_mag(0,0) = m * inv_loop_time;
    _R_mag(1,1) = m * inv_loop_time;
    _R_mag(2,2) = m * inv_loop_time;
    return;
}

//====================================Private Members==========================================//
//---------------------------------------Initialize--------------------------------------------//
void AHRS_Attitude_EKF::initialize(float g_val, float loop_time_pred_s, float loop_time_update_gyro_s, float loop_time_update_acc_s, float loop_time_update_mag_s, MatrixXf inertia_m){
    _loop_time_pred_s = loop_time_pred_s;
    _loop_time_update_gyro_s = loop_time_update_gyro_s;
    _loop_time_update_acc_s = loop_time_update_acc_s;
    _loop_time_update_mag_s = loop_time_update_mag_s;
    _g_val = g_val;
    _P.identity();
    _P *= 200.f;
    _J = inertia_m;
    _Ji = !_J;
    u(8,0) = - _g_val;
    u(9,0) = 1.f;
    _Q(0,0) = AHRS_ATTITUDE_EKF_R_S_NOISE * _loop_time_pred_s;
    _Q(1,1) = AHRS_ATTITUDE_EKF_R_S_NOISE * _loop_time_pred_s;
    _Q(2,2) = AHRS_ATTITUDE_EKF_R_S_NOISE * _loop_time_pred_s;
    _Q(3,3) = AHRS_ATTITUDE_EKF_R_A_NOISE * _loop_time_pred_s;
    _Q(4,4) = AHRS_ATTITUDE_EKF_R_A_NOISE * _loop_time_pred_s;
    _Q(5,5) = AHRS_ATTITUDE_EKF_R_A_NOISE * _loop_time_pred_s;
    _Q(6,6) = AHRS_ATTITUDE_EKF_A_NOISE * _loop_time_pred_s;
    _Q(7,7) = AHRS_ATTITUDE_EKF_A_NOISE * _loop_time_pred_s;
    _Q(8,8) = AHRS_ATTITUDE_EKF_A_NOISE * _loop_time_pred_s;
    _Q(9,9) = AHRS_ATTITUDE_EKF_M_NOISE * _loop_time_pred_s;
    _Q(10,10) = AHRS_ATTITUDE_EKF_M_NOISE * _loop_time_pred_s;
    _Q(11,11) = AHRS_ATTITUDE_EKF_M_NOISE * _loop_time_pred_s;
    _R_gyro(0,0) = AHRS_ATTITUDE_EKF_GYRO_NOISE / _loop_time_update_gyro_s;
    _R_gyro(1,1) = AHRS_ATTITUDE_EKF_GYRO_NOISE / _loop_time_update_gyro_s;
    _R_gyro(2,2) = AHRS_ATTITUDE_EKF_GYRO_NOISE / _loop_time_update_gyro_s;
    _R_accel(0,0) = AHRS_ATTITUDE_EKF_ACCEL_NOISE / _loop_time_update_acc_s;
    _R_accel(1,1) = AHRS_ATTITUDE_EKF_ACCEL_NOISE / _loop_time_update_acc_s;
    _R_accel(2,2) = AHRS_ATTITUDE_EKF_ACCEL_NOISE / _loop_time_update_acc_s;
    _R_mag(0,0) = AHRS_ATTITUDE_EKF_MAG_NOISE / _loop_time_update_mag_s;
    _R_mag(1,1) = AHRS_ATTITUDE_EKF_MAG_NOISE / _loop_time_update_mag_s;
    _R_mag(2,2) = AHRS_ATTITUDE_EKF_MAG_NOISE / _loop_time_update_mag_s;
    return;
    
}


