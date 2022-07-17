//  IMU_EKF_AV.cpp
//
//
//  Created by Andrea Vivani on 21/4/15.
//  Copyright (c) 2015 Andrea Vivani. All rights reserved.
//
#include "IMU_EKF_AV.h"
//=====================================Constructors==========================================//
IMU_EKF_AV::IMU_EKF_AV(float g_val, float loop_time_s):
u(7,1),_A(7,7),_B(7,6),_C(2,7),_P(7,7),_W(6,6),_R(2,2),_M(2,2),_K(7,2){
    _loop_time_s = loop_time_s;
    _g = g_val;
    //initialize matrices
    u(0,0) = IMU_EKF_AV_PHI0;
    u(1,0) = IMU_EKF_AV_THETA0;
    u(5,0) = IMU_EKF_AV_C_DAMP0;
    _W(0,0) = IMU_EKF_AV_GXY_NOISE * _loop_time_s;
    _W(1,1) = IMU_EKF_AV_GXY_NOISE * _loop_time_s;
    _W(2,2) = IMU_EKF_AV_GZ_NOISE * _loop_time_s;
    _W(3,3) = IMU_EKF_AV_AZ_NOISE * _loop_time_s;
    _W(4,4) = IMU_EKF_AV_C_DAMP_NOISE * _loop_time_s;
    _W(5,5) = IMU_EKF_AV_B_AZ_NOISE * _loop_time_s;
    _R(0,0) = IMU_EKF_AV_AXY_NOISE / _loop_time_s;
    _R(1,1) = IMU_EKF_AV_AXY_NOISE / _loop_time_s;
    _r_vxy = IMU_EKF_AV_VXY_NOISE;
    _r_vz = IMU_EKF_AV_VZ_NOISE;
    _r_vd = IMU_EKF_AV_VD_NOISE;
    
    //Set Angles
    Roll=u(0,0);
    Pitch=u(1,0);
    return;
}

//====================================Public Members==========================================//
//-------------------------------------g_xediction---------------------------------------------//
void IMU_EKF_AV::prediction(float g_x, float g_y, float g_z, float a_z){
    
    //Trig functions
    float sphi = sinf(u(0,0));
    float cphi = cosf(u(0,0));
    float stheta = sinf(u(1,0));
    float ctheta = cosf(u(1,0));
    float inv_ctheta = 1.0f / ctheta;
    float ttheta = stheta * inv_ctheta;
    float tmp1 = sphi * g_y + cphi * g_z;
    float tmp2 = cphi * g_y - sphi * g_z;
    
    //A matrix
    //_A.zeros(); //zeros or not?
    _A(0,0) = 1.0f + _loop_time_s * tmp2 * ttheta;
    _A(0,1) = _loop_time_s * tmp1 * inv_ctheta * inv_ctheta;
    _A(1,0) = -_loop_time_s * tmp1;
    _A(1,1) = 1.0f;
    _A(2,1) = -_loop_time_s * _g * ctheta;
    _A(2,2) = 1.0f-_loop_time_s * u(5,0);
    _A(2,3) = _loop_time_s * g_z;
    _A(2,4) = -_loop_time_s * g_y;
    _A(2,5) = -_loop_time_s * u(2,0);
    _A(3,0) = _loop_time_s * _g * cphi * ctheta;
    _A(3,1) = -_loop_time_s * _g * sphi * stheta;
    _A(3,2) = -_loop_time_s * g_z;
    _A(3,3) = 1.0f-_loop_time_s * u(5,0);
    _A(3,4) = _loop_time_s * g_x;
    _A(3,5) = -_loop_time_s * u(3,0);
    _A(4,0) = _loop_time_s * (u(6,0) - _g) * ctheta * sphi;
    _A(4,1) = _loop_time_s * (u(6,0) - _g) * cphi * stheta;
    _A(4,4) = 1.0f;
    _A(4,6) = -_loop_time_s * cphi * ctheta;
    _A(5,5) = 1.0f;
    _A(6,6) = 1.0f;
    
    //B matrix
    //_B.zeros(); //zeros or not?
    _B(0,0) = 1.0f;
    _B(0,1) = sphi * ttheta;
    _B(0,2) = cphi * ttheta;
    _B(1,1) = cphi;
    _B(1,2) = -sphi;
    _B(2,1) = -u(4,0);
    _B(2,2) = u(3,0);
    _B(3,0) = u(4,0);
    _B(3,2) = -u(2,0);
    _B(4,3) = 1.0f;
    _B(5,4) = 1.0f;
    _B(6,5) = 1.0f;
    
    //Predicted P matrix
    //P_m = A * P_p * (~A) + B * W * (~B);
    //Q=A*B*W*(~B)*(~A); //with continuous-time A (Ad=I+A*dt), it should be Q=A*B*W*(~B)*(~A)*T_samp but T_samp is already included in _W
    //Q=B*W*(~B); //it should be Q=B*W*(~B)*T_samp but T_samp is already included in _W
    _P = QuadProd(_A, _P) + QuadProd(_B, _W);
    
    //Predict state
    u(0,0) += _loop_time_s * (g_x + tmp1 * ttheta);
    u(1,0) += _loop_time_s * tmp2;
    float delta_u2 = _loop_time_s * (u(3,0) * g_z - u(4,0) * g_y - u(5,0) * u(2,0) - _g * stheta);
    float delta_u3 = _loop_time_s * (u(4,0) * g_x - u(2,0) * g_z - u(5,0) * u(3,0) + _g * sphi * ctheta);
    u(2,0) += delta_u2;
    u(3,0) += delta_u3;
    u(4,0) += _loop_time_s * (a_z + (_g - u(6,0)) * cphi * ctheta);
    /*u(5,0) += 0;
     u(6,0) += 0;*/
    
    //Set Angles
    Roll=u(0,0);
    Pitch=u(1,0);
    
    return;
}

//---------------------------------Update with accel----------------------------------------//
void IMU_EKF_AV::update(float g_x, float g_y, float g_z, float a_x, float a_y){
    
    //C matrix
    //_C.zeros(); //zeros or not?
    _C(0,2) = -u(5,0);
    _C(0,3) = g_z;
    _C(0,4) = -g_y;
    _C(0,5) = -u(2,0);
    _C(1,2) = -g_z;
    _C(1,3) = -u(5,0);
    _C(1,4) = g_x;
    _C(1,5) = -u(3,0);
    
    //Delta measures
    MatrixXf Delta_m(2,1);
    Delta_m(0,0) = a_x + u(4,0) * g_y - u(3,0) * g_z + u(5,0) * u(2,0);
    Delta_m(1,0) = a_y + u(2,0) * g_z - u(4,0) * g_x + u(5,0) * u(3,0);
    
    //Gain matrix K
    _M = QuadProd(_C,_P) + _R;
    _K = _P * (~_C) * (!_M);
    
    //Correct state vector
    u += _K * Delta_m;
    
    //Updated P matrix
    _P -= _K * _C * _P;
    //_P=(_P+(~_P))*0.5; //guarantees P to be symmetric
    
    //Set Angles
    Roll=u(0,0);
    Pitch=u(1,0);

    return;
}

//-------------------------Update with velocity along x,y local------------------------------//
void IMU_EKF_AV::update_vel_xy(float v_x, float v_y, float dt_s){
    
    //R matrix
    MatrixXf R_tmp(2,2);
    R_tmp(0,0) = _r_vxy / dt_s;
    R_tmp (1,1) = R_tmp(0,0);
    
    //C matrix
    MatrixXf C_tmp(2,7);
    //C_tmp.zeros(); //zeros or not?
    C_tmp(0,2) = 1.f;
    C_tmp(1,3) = 1.f;
    
    //Delta measures
    MatrixXf Delta_m(2,1);
    Delta_m(0,0) = v_x - u(2,0);
    Delta_m(1,0) = v_y - u(3,0);
    
    //Gain matrix K
    _M = QuadProd(C_tmp,_P) + R_tmp;
    _K = _P * (~C_tmp) * (!_M);
    
    //Correct state vector
    u += _K * Delta_m;
    
    //Updated P matrix
    _P -= _K * C_tmp * _P;
    //_P=(_P+(~_P))*0.5; //guarantees P to be symmetric
    
    //Set Angles
    Roll=u(0,0);
    Pitch=u(1,0);

    return;
}

//--------------------------Update with velocity along z local--------------------------------//
void IMU_EKF_AV::update_vel_z(float v_z, float dt_s){
    
    
    //C matrix
    MatrixXf C_tmp(1,7);
    //C_tmp.zeros(); //zeros or not?
    C_tmp(0,4) = 1.f;
    
    //Delta measures
    float Delta_m = v_z - u(4,0);

    //Gain matrix K (change K with _K in correction and update of P)
    /*_M = QuadProd(C_tmp,_P) + (_r_vz / dt_s);
    _K = _P * (~C_tmp) * (!_M);*/
    
     
     //Faster Gain matrix K
    float inv_m = 1.f / (_P(4,4) + (_r_vz / dt_s));
    MatrixXf K(7,1);
    K(0,0) = _P(0,4) * inv_m;
    K(1,0) = _P(1,4) * inv_m;
    K(2,0) = _P(2,4) * inv_m;
    K(3,0) = _P(3,4) * inv_m;
    K(4,0) = _P(4,4) * inv_m;
    K(5,0) = _P(5,4) * inv_m;
    K(6,0) = _P(6,4) * inv_m;

    //Correct state vector
    u += K * Delta_m;
    
    //Updated P matrix
    _P -= K * C_tmp * _P;
    //_P=(_P+(~_P))*0.5; //guarantees P to be symmetric
    
    
    //Set Angles
    Roll=u(0,0);
    Pitch=u(1,0);
    
    return;
}

//--------------------------Update with velocity along d global--------------------------------//
void IMU_EKF_AV::update_vel_d(float v_d, float dt_s){
    
    //Trig functions
    float sphi = sinf(u(0,0));
    float cphi = cosf(u(0,0));
    float stheta = sinf(u(1,0));
    float ctheta = cosf(u(1,0));
    
    //C matrix
    MatrixXf C_tmp(1,7);
    //C_tmp.zeros(); //zeros or not?
    C_tmp(0,0) = u(3,0) * cphi * ctheta - u(4,0) * ctheta * sphi;
    C_tmp(0,1) = - u(2,0) * ctheta - u(4,0) * cphi * stheta - u(3,0) * sphi * stheta;
    C_tmp(0,2) = - stheta;
    C_tmp(0,3) = ctheta * sphi;
    C_tmp(0,4) = cphi * ctheta;
    
    //Delta measures
    float Delta_m = v_d - (u(4,0) * cphi * ctheta - u(2,0) * stheta + u(3,0) * ctheta * sphi);
    
    //Gain matrix K
    _M = QuadProd(C_tmp,_P) + (_r_vd / dt_s);
    _K = _P * (~C_tmp) * (!_M);
    
    //Correct state vector
    u += _K * Delta_m;
    
    //Updated P matrix
    _P -= _K * C_tmp * _P;
    //_P=(_P+(~_P))*0.5; //guarantees P to be symmetric
    
    //Set Angles
    Roll=u(0,0);
    Pitch=u(1,0);
    
    return;
}

//-----------------------------------Starting values-----------------------------------------//
void IMU_EKF_AV::set_starting_values(float Phi_0, float Theta_0, float c_damp_0){
    u(0,0) = Phi_0;
    u(1,0) = Theta_0;
    u(5,0) = c_damp_0;
    
    //Set Angles
    Roll=u(0,0);
    Pitch=u(1,0);
    
    return;
}

//-------------------------------------Input noises------------------------------------------//
void IMU_EKF_AV::set_input_noises(float g_xy, float g_z, float a_z, float c_damp, float b_az){
    _W(0,0) = g_xy * _loop_time_s;
    _W(1,1) = g_xy * _loop_time_s;
    _W(2,2) = g_z * _loop_time_s;
    _W(3,3) = a_z * _loop_time_s;
    _W(4,4) = c_damp * _loop_time_s;
    _W(5,5) = b_az * _loop_time_s;
    
    return;
}

//----------------------------------Output accel noises--------------------------------------//
void IMU_EKF_AV::set_acc_noise(float a){
    float inv_loop_time_s = 1.0f / _loop_time_s;
    _R(0,0) = a * inv_loop_time_s;
    _R(1,1) = a * inv_loop_time_s;
    
    return;
}
