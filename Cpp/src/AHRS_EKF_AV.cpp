//  AHRS_EKF_AV.cpp
//
//
//  Created by Andrea Vivani on 21/4/15.
//  Copyright (c) 2015 Andrea Vivani. All rights reserved.
//

#include "AHRS_EKF_AV.h"

//=====================================Constructors==========================================//
AHRS_EKF_AV::AHRS_EKF_AV(float g_val, float loop_time_pred_update_acc_s) : AHRS_EKF_AV(g_val, loop_time_pred_update_acc_s, loop_time_pred_update_acc_s){};

AHRS_EKF_AV::AHRS_EKF_AV(float g_val, float loop_time_pred_update_acc_s, float loop_time_mag_s):
u(12,1),_A(12,12),_B(12,10),_C_acc(2,12),_C_mag(3,12),_P(12,12),_W(10,10),_R_acc(2,2),_R_mag(3,3),_M(1,1),_K(1,1){
    _loop_time_pred_update_acc_s = loop_time_pred_update_acc_s;
    _loop_time_mag_s = loop_time_mag_s;
    _g = g_val;
    //initialize matrices
    u(0,0) = AHRS_EKF_AV_PHI0;
    u(1,0) = AHRS_EKF_AV_THETA0;
    u(2,0) = AHRS_EKF_AV_PSI0;
    u(6,0) = AHRS_EKF_AV_C_DAMP0;
    u(8,0) = AHRS_EKF_AV_INCL0;
    _W(0,0) = AHRS_EKF_AV_GXY_NOISE * _loop_time_pred_update_acc_s;
    _W(1,1) = AHRS_EKF_AV_GXY_NOISE * _loop_time_pred_update_acc_s;
    _W(2,2) = AHRS_EKF_AV_GZ_NOISE * _loop_time_pred_update_acc_s;
    _W(3,3) = AHRS_EKF_AV_AZ_NOISE * _loop_time_pred_update_acc_s;
    _W(4,4) = AHRS_EKF_AV_C_DAMP_NOISE * _loop_time_pred_update_acc_s;
    _W(5,5) = AHRS_EKF_AV_B_AZ_NOISE * _loop_time_pred_update_acc_s;
    _W(6,6) = AHRS_EKF_AV_INCL_NOISE * _loop_time_pred_update_acc_s;
    _W(7,7) = AHRS_EKF_AV_B_G_NOISE * _loop_time_pred_update_acc_s;
    _W(8,8) = AHRS_EKF_AV_B_G_NOISE * _loop_time_pred_update_acc_s;
    _W(9,9) = AHRS_EKF_AV_B_G_NOISE * _loop_time_pred_update_acc_s;
    _R_acc(0,0) = AHRS_EKF_AV_A_NOISE / _loop_time_pred_update_acc_s;
    _R_acc(1,1) = AHRS_EKF_AV_A_NOISE / _loop_time_pred_update_acc_s;
    _R_mag(0,0) = AHRS_EKF_AV_M_NOISE / _loop_time_mag_s;
    _R_mag(1,1) = AHRS_EKF_AV_M_NOISE / _loop_time_mag_s;
    _R_mag(2,2) = AHRS_EKF_AV_M_NOISE / _loop_time_mag_s;
    _r_vxy = AHRS_EKF_AV_VXY_NOISE;
    _r_vz = AHRS_EKF_AV_VZ_NOISE;
    _r_vne = AHRS_EKF_AV_VNE_NOISE;
    _r_vd = AHRS_EKF_AV_VD_NOISE;
    
    //Set Angles
    Roll = u(0,0);
    Pitch = u(1,0);
    Yaw = u(2,0);
    
    return;
}

//====================================Public Members==========================================//
//-------------------------------------Prediction---------------------------------------------//
void AHRS_EKF_AV::prediction(float g_x, float g_y, float g_z, float a_z){
    
    //Trig functions
    float sphi = sinf(u(0,0));
    float cphi = cosf(u(0,0));
    float stheta = sinf(u(1,0));
    float ctheta = cosf(u(1,0));
    float inv_ctheta = 1.0f / ctheta;
    float ttheta = stheta * inv_ctheta;
    float pr = g_x - u(9,0);
    float qr = g_y - u(10,0);
    float rr = g_z - u(11,0);
    float tmp1 = sphi * qr + cphi * rr;
    float tmp2 = cphi * qr - sphi * rr;
    
    //A matrix
    //_A.zeros(); //zeros or not?
    _A(0,0) = 1.0f + _loop_time_pred_update_acc_s * tmp2 * ttheta;
    _A(0,1) = _loop_time_pred_update_acc_s * tmp1 * inv_ctheta * inv_ctheta;
    _A(0,9) = - _loop_time_pred_update_acc_s;
    _A(0,10) = - _loop_time_pred_update_acc_s * sphi * ttheta;
    _A(0,11) = - _loop_time_pred_update_acc_s * cphi * ttheta;
    _A(1,0) = - _loop_time_pred_update_acc_s * tmp1;
    _A(1,1) = 1.0f;
    _A(1,10) = - _loop_time_pred_update_acc_s * cphi;
    _A(1,11) = _loop_time_pred_update_acc_s * sphi;
    _A(2,0) = _loop_time_pred_update_acc_s * tmp2 * inv_ctheta;
    _A(2,1) = _loop_time_pred_update_acc_s * tmp1 * ttheta * inv_ctheta;
    _A(2,2) = 1.0f;
    _A(2,10) = - _loop_time_pred_update_acc_s * sphi * inv_ctheta;
    _A(2,11) = - _loop_time_pred_update_acc_s * cphi * inv_ctheta;
    _A(3,1) = - _loop_time_pred_update_acc_s * _g * ctheta;
    _A(3,3) = 1.0f - _loop_time_pred_update_acc_s * u(6,0);
    _A(3,4) = _loop_time_pred_update_acc_s * rr;
    _A(3,5) = - _loop_time_pred_update_acc_s * qr;
    _A(3,6) = - _loop_time_pred_update_acc_s * u(3,0);
    _A(4,0) = _loop_time_pred_update_acc_s * _g * cphi * ctheta;
    _A(4,1) = - _loop_time_pred_update_acc_s * _g * sphi * stheta;
    _A(4,3) = - _loop_time_pred_update_acc_s * rr;
    _A(4,4) = 1.0f-_loop_time_pred_update_acc_s * u(6,0);
    _A(4,5) = _loop_time_pred_update_acc_s * pr;
    _A(4,6) = - _loop_time_pred_update_acc_s * u(4,0);
    _A(5,0) = _loop_time_pred_update_acc_s * (u(7,0) - _g) * ctheta * sphi;
    _A(5,1) = _loop_time_pred_update_acc_s * (u(7,0) - _g) * cphi * stheta;
    _A(5,5) = 1.0f;
    _A(5,7) = - _loop_time_pred_update_acc_s * cphi * ctheta;
    _A(6,6) = 1.0f;
    _A(7,7) = 1.0f;
    _A(8,8) = 1.0f;
    _A(9,9) = 1.0f;
    _A(10,10) = 1.0f;
    _A(11,11) = 1.0f;
    
    
    //B matrix
    //_B.zeros(); //zeros or not?
    _B(0,0) = 1.0f;
    _B(0,1) = sphi * ttheta;
    _B(0,2) = cphi * ttheta;
    _B(1,1) = cphi;
    _B(1,2) = -sphi;
    _B(2,1) = sphi * inv_ctheta;
    _B(2,2) = cphi * inv_ctheta;
    _B(3,1) = -u(5,0);
    _B(3,2) = u(4,0);
    _B(4,0) = u(5,0);
    _B(4,2) = -u(3,0);
    _B(5,3) = 1.0f;
    _B(6,4) = 1.0f;
    _B(7,5) = 1.0f;
    _B(8,6) = 1.0f;
    _B(9,7) = 1.0f;
    _B(10,8) = 1.0f;
    _B(11,9) = 1.0f;
    
    //Predicted P matrix
    //P_m = A * P_p * (~A) + B * W * (~B);
    //Q=A*B*W*(~B)*(~A); //with continuous-time A (Ad=I+A*dt), it should be Q=A*B*W*(~B)*(~A)*T_samp but T_samp is already included in _W
    //Q=B*W*(~B); //it should be Q=B*W*(~B)*T_samp but T_samp is already included in _W
    _P = QuadProd(_A, _P) + QuadProd(_B, _W);
    
    //Predict state
    u(0,0) += _loop_time_pred_update_acc_s * (pr + tmp1 * ttheta);
    u(1,0) += _loop_time_pred_update_acc_s * tmp2;
    u(2,0) += _loop_time_pred_update_acc_s * tmp1 * inv_ctheta;
    float delta_u3 = _loop_time_pred_update_acc_s * (u(4,0) * rr - u(5,0) * qr - u(6,0) * u(3,0) - _g * stheta);
    float delta_u4 = _loop_time_pred_update_acc_s * (u(5,0) * pr - u(3,0) * rr - u(6,0) * u(4,0) + _g * sphi * ctheta);
    u(3,0) += delta_u3;
    u(4,0) += delta_u4;
    u(5,0) += _loop_time_pred_update_acc_s * (a_z + (_g - u(7,0)) * cphi * ctheta);
    /*u(6,0) += 0;
     u(7,0) += 0;
     u(8,0) += 0;
     u(9,0) += 0;
     u(10,0) += 0;
     u(11,0) += 0;*/
    
    //Set Angles
    Roll=u(0,0);
    Pitch=u(1,0);
    Yaw=u(2,0);
    
    return;
}

//---------------------------------Update with accel----------------------------------------//
void AHRS_EKF_AV::update_accel(float g_x, float g_y, float g_z, float a_x, float a_y){
    
    //Remove bias
    float pr = g_x - u(9,0);
    float qr = g_y - u(10,0);
    float rr = g_z - u(11,0);
    
    //C matrix
    //_C_acc.zeros(); //zeros or not?
    _C_acc(0,3) = - u(6,0);
    _C_acc(0,4) = rr;
    _C_acc(0,5) = - qr;
    _C_acc(0,6) = - u(3,0);
    _C_acc(1,3) = - rr;
    _C_acc(1,4) = - u(6,0);
    _C_acc(1,5) = pr;
    _C_acc(1,6) = - u(4,0);

    //Delta measures
    MatrixXf Delta_m(2,1);
    Delta_m(0,0) = a_x + u(5,0) * qr - u(4,0) * rr + u(6,0) * u(3,0);
    Delta_m(1,0) = a_y + u(3,0) * rr - u(5,0) * pr + u(6,0) * u(4,0);
    
    //Gain matrix K
    _M = QuadProd(_C_acc,_P) + _R_acc;
    _K = _P * (~_C_acc) * (!_M);

    //Correct state vector
    u += _K * Delta_m;
    
    //Updated P matrix
    _P -= _K * _C_acc * _P;
    //_P=(_P+(~_P))*0.5; //guarantees P to be symmetric
    
    //Set Angles
    Roll = u(0,0);
    Pitch = u(1,0);
    Yaw = u(2,0);
    
    return;
}

//-----------------------------------Update with mag-----------------------------------------//
void AHRS_EKF_AV::update_mag(float m_x, float m_y, float m_z){
    
    //Normalize readings
    float inv_norm=1.f / sqrtf(m_x * m_x + m_y * m_y + m_z * m_z);
    if (isnan(inv_norm)||isinf(inv_norm)){
        inv_norm=1.f;
    }
    m_x *= inv_norm;
    m_y *= inv_norm;
    m_z *= inv_norm;
    
    //Trig functions
    float sinc = sinf(u(8,0));
    float cinc = cosf(u(8,0));
    float sphi = sinf(u(0,0));
    float cphi = cosf(u(0,0));
    float stheta = sinf(u(1,0));
    float ctheta = cosf(u(1,0));
    float spsi = sinf(u(2,0));
    float cpsi = cosf(u(2,0));
    
    //C matrix
    //_C_mag.zeros(); //zeros or not?
    float tmp1 = cphi * spsi - cpsi * sphi * stheta;
    _C_mag(0,1) =- sinc * ctheta - cinc * cpsi * stheta;
    _C_mag(0,2) = - cinc * ctheta * spsi;
    _C_mag(0,7) = - cinc * stheta - cpsi * ctheta * sinc;
    _C_mag(1,0) = cinc * (sphi * spsi + cphi * cpsi * stheta) + cphi * sinc * ctheta;
    _C_mag(1,1) = cinc * cpsi * ctheta * sphi - sinc * sphi * stheta;
    _C_mag(1,2) = - cinc * (cphi * cpsi + sphi * spsi * stheta);
    _C_mag(1,8) = sinc * tmp1 + cinc * ctheta * sphi;
    _C_mag(2,0) = cinc * tmp1 - sinc * ctheta * sphi;
    _C_mag(2,1) = cinc * cphi * cpsi * ctheta - cphi * sinc * stheta;
    _C_mag(2,2) = cinc * (cpsi * sphi - cphi * spsi * stheta);
    _C_mag(2,8) = cinc * cphi * ctheta - sinc * (sphi * spsi + cphi * cpsi * stheta);

    MatrixXf Delta_m(3,1);
    Delta_m(0,0) = m_x - (cinc * cpsi * ctheta - sinc * stheta);
    Delta_m(1,0) = m_y - (sinc * ctheta * sphi - cinc * tmp1);
    Delta_m(2,0) = m_z - (cinc * (sphi * spsi + cphi * cpsi * stheta) + cphi * sinc * ctheta);
    
    //Gain matrix K
    _M = QuadProd(_C_mag,_P) + _R_mag;
    _K = _P * (~_C_mag) * (!_M);
    
    //Correct state vector
    u += _K * Delta_m;
    
    //Updated P matrix
    _P -= _K * _C_mag * _P;
    //_P=(_P+(~_P))*0.5; //guarantees P to be symmetric
    
    //Set Angles
    Roll = u(0,0);
    Pitch = u(1,0);
    Yaw = u(2,0);
    
    return;
}

//-------------------------Update with velocity along x,y local------------------------------//
void AHRS_EKF_AV::update_vel_xy(float v_x, float v_y, float dt_s){
    
    //R matrix
    MatrixXf R_tmp(2,2);
    R_tmp(0,0) = _r_vxy / dt_s;
    R_tmp (1,1) = R_tmp(0,0);
    
    //C matrix
    MatrixXf C_tmp(2,12);
    //C_tmp.zeros(); //zeros or not?
    C_tmp(0,3) = 1.f;
    C_tmp(1,4) = 1.f;
    
    //Delta measures
    MatrixXf Delta_m(2,1);
    Delta_m(0,0) = v_x - u(3,0);
    Delta_m(1,0) = v_y - u(4,0);
    
    //Gain matrix K
    _M = QuadProd(C_tmp,_P) + R_tmp;
    _K = _P * (~C_tmp) * (!_M);
    
    //Correct state vector
    u += _K * Delta_m;
    
    //Updated P matrix
    _P -= _K * C_tmp * _P;
    //_P=(_P+(~_P))*0.5; //guarantees P to be symmetric
    
    //Set Angles
    Roll = u(0,0);
    Pitch = u(1,0);
    
    return;
}

//--------------------------Update with velocity along z local--------------------------------//
void AHRS_EKF_AV::update_vel_z(float v_z, float dt_s){
    
    
    //C matrix
    MatrixXf C_tmp(1,12);
    //C_tmp.zeros(); //zeros or not?
    C_tmp(0,5) = 1.f;
    
    //Delta measures
    float Delta_m = v_z - u(5,0);
    
    //Gain matrix K (change K with _K in correction and update of P)
    /*_M = QuadProd(C_tmp,_P) + (_r_vz / dt_s);
     _K = _P * (~C_tmp) * (!_M);*/
    
    
    //Faster Gain matrix K
    float inv_m = 1.f / (_P(5,5) + (_r_vz / dt_s));
    MatrixXf K(12,1);
    K(0,0) = _P(0,5) * inv_m;
    K(1,0) = _P(1,5) * inv_m;
    K(2,0) = _P(2,5) * inv_m;
    K(3,0) = _P(3,5) * inv_m;
    K(4,0) = _P(4,5) * inv_m;
    K(5,0) = _P(5,5) * inv_m;
    K(6,0) = _P(6,5) * inv_m;
    K(7,0) = _P(7,5) * inv_m;
    K(8,0) = _P(8,5) * inv_m;
    K(9,0) = _P(9,5) * inv_m;
    K(10,0) = _P(10,5) * inv_m;
    K(11,0) = _P(11,5) * inv_m;
    
    //Correct state vector
    u += K * Delta_m;
    
    //Updated P matrix
    _P -= K * C_tmp * _P;
    //_P=(_P+(~_P))*0.5; //guarantees P to be symmetric
    
    
    //Set Angles
    Roll = u(0,0);
    Pitch = u(1,0);
    
    return;
}

//-------------------------Update with velocity along n,e global------------------------------//
void AHRS_EKF_AV::update_vel_ne(float v_n, float v_e, float dt_s){
    
    //R matrix
    MatrixXf R_tmp(2,2);
    R_tmp(0,0) = _r_vne / dt_s;
    R_tmp (1,1) = R_tmp(0,0);
    
    //Trig functions
    float sphi = sinf(u(0,0));
    float cphi = cosf(u(0,0));
    float stheta = sinf(u(1,0));
    float ctheta = cosf(u(1,0));
    float spsi = sinf(u(2,0));
    float cpsi = cosf(u(2,0));
    
    //C matrix
    MatrixXf C_tmp(2,12);
    //C_tmp.zeros(); //zeros or not?
    float tmp1 = sphi * spsi + cphi * cpsi * stheta;
    float tmp2 = cphi * cpsi + sphi * spsi * stheta;
    float tmp3 = cphi * spsi - cpsi * sphi * stheta;
    float tmp4 = cpsi * sphi - cphi * spsi * stheta;
    C_tmp(0,0) = u(4,0) * tmp1 + u(5,0) * tmp3;
    C_tmp(0,1) = u(5,0) * cphi * cpsi * ctheta - u(3,0) * cpsi * stheta + u(4,0) * cpsi * ctheta * sphi;
    C_tmp(0,2) = u(5,0) *  - u(4,0)* tmp2 - u(3,0) * ctheta * spsi;
    C_tmp(0,3) = cpsi * ctheta;
    C_tmp(0,4) = - tmp3;
    C_tmp(0,5) = tmp1;
    C_tmp(1,0) = - u(4,0) * tmp4 - u(5,0) * tmp2;
    C_tmp(1,1) = u(5,0) * cphi * ctheta * spsi - u(3,0) * spsi * stheta + u(4,0) * ctheta * sphi * spsi;
    C_tmp(1,2) = u(5,0) * tmp1 - u(4,0) * tmp3 + u(3,0) * cpsi * ctheta;
    C_tmp(1,3) = ctheta * spsi;
    C_tmp(1,4) = tmp2;
    C_tmp(1,5) = tmp4;
    
    //Delta measures
    MatrixXf Delta_m(2,1);
    Delta_m(0,0) = v_n - (u(5,0) * tmp1 - u(4,0) * tmp3 + u(3,0) * cpsi * ctheta);
    Delta_m(1,0) = v_e - (u(4,0) * tmp2 - u(5,0) * tmp4 + u(3,0) * ctheta * spsi);
    
    //Gain matrix K
    _M = QuadProd(C_tmp,_P) + R_tmp;
    _K = _P * (~C_tmp) * (!_M);
    
    //Correct state vector
    u += _K * Delta_m;
    
    //Updated P matrix
    _P -= _K * C_tmp * _P;
    //_P=(_P+(~_P))*0.5; //guarantees P to be symmetric
    
    //Set Angles
    Roll = u(0,0);
    Pitch = u(1,0);
    
    return;
}

//--------------------------Update with velocity along d global--------------------------------//
void AHRS_EKF_AV::update_vel_d(float v_d, float dt_s){
    
    //Trig functions
    float sphi = sinf(u(0,0));
    float cphi = cosf(u(0,0));
    float stheta = sinf(u(1,0));
    float ctheta = cosf(u(1,0));
    
    //C matrix
    MatrixXf C_tmp(1,12);
    //C_tmp.zeros(); //zeros or not?
    C_tmp(0,0) = u(4,0) * cphi * ctheta - u(5,0) * ctheta * sphi;
    C_tmp(0,1) = - u(3,0) * ctheta - u(5,0) * cphi * stheta - u(4,0) * sphi * stheta;
    C_tmp(0,3) = - stheta;
    C_tmp(0,4) = ctheta * sphi;
    C_tmp(0,5) = cphi * ctheta;
    
    
    //Delta measures
    float Delta_m = v_d - (u(5,0) * cphi * ctheta - u(3,0) * stheta + u(4,0) * ctheta * sphi);
    
    //Gain matrix K
    _M = QuadProd(C_tmp,_P) + (_r_vd / dt_s);
    _K = _P * (~C_tmp) * (!_M);
    
    //Correct state vector
    u += _K * Delta_m;
    
    //Updated P matrix
    _P -= _K * C_tmp * _P;
    //_P=(_P+(~_P))*0.5; //guarantees P to be symmetric
    
    //Set Angles
    Roll = u(0,0);
    Pitch = u(1,0);
    
    return;
}

//-----------------------------------Starting values-----------------------------------------//
void AHRS_EKF_AV::set_starting_values(float Phi_0, float Theta_0, float Psi_0, float c_damp_0, float incl_0){
    u(0,0) = Phi_0;
    u(1,0) = Theta_0;
    u(2,0) = Psi_0;
    u(6,0) = c_damp_0;
    u(8,0) = incl_0;
    
    //Set Angles
    Roll = u(0,0);
    Pitch = u(1,0);
    Yaw = u(2,0);
    
    return;
}

//-------------------------------------Input noises------------------------------------------//
void AHRS_EKF_AV::set_input_noises(float g_xy, float g_z, float a_z, float c_damp, float b_az, float incl, float b_g){
    _W(0,0) = g_xy * _loop_time_pred_update_acc_s;
    _W(1,1) = g_xy * _loop_time_pred_update_acc_s;
    _W(2,2) = g_z * _loop_time_pred_update_acc_s;
    _W(3,3) = a_z * _loop_time_pred_update_acc_s;
    _W(4,4) = c_damp * _loop_time_pred_update_acc_s;
    _W(5,5) = b_az * _loop_time_pred_update_acc_s;
    _W(6,6) = incl * _loop_time_pred_update_acc_s;
    _W(7,7) = b_g * _loop_time_pred_update_acc_s;
    _W(8,8) = b_g * _loop_time_pred_update_acc_s;
    _W(9,9) = b_g * _loop_time_pred_update_acc_s;
    
    return;
}

//----------------------------------Output accel noises--------------------------------------//
void AHRS_EKF_AV::set_acc_noises(float a){
    float inv_loop_time_pred_update_acc_s = 1.0f / _loop_time_pred_update_acc_s;
    _R_acc(0,0) = a * inv_loop_time_pred_update_acc_s;
    _R_acc(1,1) = a * inv_loop_time_pred_update_acc_s;
    
    return;
}

//-----------------------------------Output mag noises---------------------------------------//	
void AHRS_EKF_AV::set_mag_noises(float m){
    float inv_loop_time_mag_s = 1.0f / _loop_time_mag_s;
    _R_mag(0,0) = m * inv_loop_time_mag_s;
    _R_mag(1,1) = m * inv_loop_time_mag_s;
    _R_mag(2,2) = m * inv_loop_time_mag_s;
    return;
}