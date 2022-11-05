//  IMU_EKF_AV.c
//
//
//  Created by Andrea Vivani on 25/06/22.
//  Copyright (c) 2022 Andrea Vivani. All rights reserved.
//
#include "IMU_EKF_AV.h"
#include "NumMethods.h"

//=====================================Global variables==========================================//

float _IMU_EKF_u_data[7];
float _A_data[7 * 7];
float _B_data[7 * 6];
float _C_data[2 * 7];
float _P_data[7 * 7];
float _W_data[6 * 6];
float _R_data[2 * 2];
float _M_data[2 * 2];
float _K_data[7 * 2];
Matrix IMU_EKF_u; // state vector (Roll=Phi, Pitch=Theta, Xd, Yd, Zd, c_damp, b_az) angles in rad, angular velocities in rad/s, velocities in m/s, c_damp in N*s/m
Matrix _A;  //state matrix
Matrix _B;  //input matrix
Matrix _C;  //output matrix
Matrix _P;  //errors expected value matrix
Matrix _W; //gyro and acc_z noises covariance matrix (g_x, g_y, g_z, a_z, c_damp, b_az)
Matrix _R;  //acc_x and acc_y noises covariance matrix
Matrix _M;  //temporary matrix
Matrix _K;  //gain matrix

// Temporary variables
float _TMP1_data[7 * 2];
float _TMP2_data[2 * 2];
float _TMP3_data[7 * 1];
float _TMP4_data[7 * 7];
float _TMP5_data[7 * 7];
Matrix TMP1, TMP2, TMP3, TMP4, TMP5;

float _loop_time_s; //loop time, in s
float _g; //value of gravitational acceleration, in m/s^2
float _r_vxy, _r_vz, _r_vd; //velocities noise covariances

//=====================================Initialization==========================================//
void IMU_EKF_init(AHRS_State_t *ahrs, float g_val, float loop_time_s)
{
  _loop_time_s = loop_time_s;
  _g = g_val;

  //initialize matrices
  IMU_EKF_u = newMatrix(7, 1, _IMU_EKF_u_data);
  _A = newMatrix(7, 7, _A_data);
  _B = newMatrix(7, 6, _B_data);
  _C = newMatrix(2, 7, _C_data);
  _P = newMatrix(7, 7, _P_data);
  _W = newMatrix(6, 6, _W_data);
  _R = newMatrix(2, 2, _R_data);
  _M = newMatrix(2, 2, _M_data);
  _K = newMatrix(7, 2, _K_data);
  TMP1 = newMatrix(7, 2, _TMP1_data);
  TMP2 = newMatrix(2, 2, _TMP2_data);
  TMP3 = newMatrix(7, 1, _TMP3_data);
  TMP4 = newMatrix(7, 7, _TMP4_data);
  TMP5 = newMatrix(7, 7, _TMP5_data);
  matZeros(IMU_EKF_u);
  matZeros(_A);
  matZeros(_B);
  matZeros(_C);
  matZeros(_P);
  matZeros(_W);
  matZeros(_R);
  matZeros(_M);
  matZeros(_K);
  matZeros(TMP1);
  matZeros(TMP2);
  matZeros(TMP3);
  matZeros(TMP4);
  matZeros(TMP5);
  ELEM(IMU_EKF_u, 0, 0) = IMU_EKF_AV_PHI0;
  ELEM(IMU_EKF_u, 1, 0) = IMU_EKF_AV_THETA0;
  ELEM(IMU_EKF_u, 5, 0) = IMU_EKF_AV_C_DAMP0;
  ELEM(_W, 0, 0) = IMU_EKF_AV_GXY_NOISE * _loop_time_s;
  ELEM(_W, 1, 1) = IMU_EKF_AV_GXY_NOISE * _loop_time_s;
  ELEM(_W, 2, 2) = IMU_EKF_AV_GZ_NOISE * _loop_time_s;
  ELEM(_W, 3, 3) = IMU_EKF_AV_AZ_NOISE * _loop_time_s;
  ELEM(_W, 4, 4) = IMU_EKF_AV_C_DAMP_NOISE * _loop_time_s;
  ELEM(_W, 5, 5) = IMU_EKF_AV_B_AZ_NOISE * _loop_time_s;
  ELEM(_R, 0, 0) = IMU_EKF_AV_AXY_NOISE / _loop_time_s;
  ELEM(_R, 1, 1) = IMU_EKF_AV_AXY_NOISE / _loop_time_s;
  _r_vxy = IMU_EKF_AV_VXY_NOISE;
  _r_vz = IMU_EKF_AV_VZ_NOISE;
  _r_vd = IMU_EKF_AV_VD_NOISE;

  //Set Angles
  //TODO Check if they're ok
  ahrs->e.thx = ELEM(IMU_EKF_u, 0, 0); //u(0,0) is roll according to IMU ref. frame
  ahrs->e.thy = ELEM(IMU_EKF_u, 1, 0); //u(1,0) is pitch according to IMU ref. frame

  return;
}

//=====================================Calculation============================================//
//-------------------------------------Prediction---------------------------------------------//
void IMU_EKF_prediction(float g_x, float g_y, float g_z, float a_z)
{
  float delta_u2, delta_u3;

  //Trig functions
  float sphi = sinf(ELEM(IMU_EKF_u, 0, 0));
  float cphi = cosf(ELEM(IMU_EKF_u, 0, 0));
  float stheta = sinf(ELEM(IMU_EKF_u, 1, 0));
  float ctheta = cosf(ELEM(IMU_EKF_u, 1, 0));
  float inv_ctheta = 1.0f / ctheta;
  float ttheta = stheta * inv_ctheta;
  float tmp1 = sphi * g_y + cphi * g_z;
  float tmp2 = cphi * g_y - sphi * g_z;

  //A matrix
  //_A.zeros(); //zeros or not?
  ELEM(_A, 0, 0) = 1.0f + _loop_time_s * tmp2 * ttheta;
  ELEM(_A, 0, 1) = _loop_time_s * tmp1 * inv_ctheta * inv_ctheta;
  ELEM(_A, 1, 0) = -_loop_time_s * tmp1;
  ELEM(_A, 1, 1) = 1.0f;
  ELEM(_A, 2, 1) = -_loop_time_s * _g * ctheta;
  ELEM(_A, 2, 2) = 1.0f - _loop_time_s * ELEM(IMU_EKF_u, 5, 0);
  ELEM(_A, 2, 3) = _loop_time_s * g_z;
  ELEM(_A, 2, 4) = -_loop_time_s * g_y;
  ELEM(_A, 2, 5) = -_loop_time_s * ELEM(IMU_EKF_u, 2, 0);
  ELEM(_A, 3, 0) = _loop_time_s * _g * cphi * ctheta;
  ELEM(_A, 3, 1) = -_loop_time_s * _g * sphi * stheta;
  ELEM(_A, 3, 2) = -_loop_time_s * g_z;
  ELEM(_A, 3, 3) = 1.0f - _loop_time_s * ELEM(IMU_EKF_u, 5, 0);
  ELEM(_A, 3, 4) = _loop_time_s * g_x;
  ELEM(_A, 3, 5) = -_loop_time_s * ELEM(IMU_EKF_u, 3, 0);
  ELEM(_A, 4, 0) = _loop_time_s * (ELEM(IMU_EKF_u, 6, 0) - _g) * ctheta * sphi;
  ELEM(_A, 4, 1) = _loop_time_s * (ELEM(IMU_EKF_u, 6, 0) - _g) * cphi * stheta;
  ELEM(_A, 4, 4) = 1.0f;
  ELEM(_A, 4, 6) = -_loop_time_s * cphi * ctheta;
  ELEM(_A, 5, 5) = 1.0f;
  ELEM(_A, 6, 6) = 1.0f;

  //B matrix
  //_B.zeros(); //zeros or not?
  ELEM(_B, 0, 0) = 1.0f;
  ELEM(_B, 0, 1) = sphi * ttheta;
  ELEM(_B, 0, 2) = cphi * ttheta;
  ELEM(_B, 1, 1) = cphi;
  ELEM(_B, 1, 2) = -sphi;
  ELEM(_B, 2, 1) = -ELEM(IMU_EKF_u, 4, 0);
  ELEM(_B, 2, 2) = ELEM(IMU_EKF_u, 3, 0);
  ELEM(_B, 3, 0) = ELEM(IMU_EKF_u, 4, 0);
  ELEM(_B, 3, 2) = -ELEM(IMU_EKF_u, 2, 0);
  ELEM(_B, 4, 3) = 1.0f;
  ELEM(_B, 5, 4) = 1.0f;
  ELEM(_B, 6, 5) = 1.0f;

  //Predicted P matrix
  //P_m = A * P_p * (~A) + B * W * (~B);
  //Q=A*B*W*(~B)*(~A); //with continuous-time A (Ad=I+A*dt), it should be Q=A*B*W*(~B)*(~A)*T_samp but T_samp is already included in _W
  //Q=B*W*(~B); //it should be Q=B*W*(~B)*T_samp but T_samp is already included in _W
  //_P = QuadProd(_A, _P) + QuadProd(_B, _W);
  QuadProd(_A, _P, TMP4);
  matCopyData(TMP4, _P);
  QuadProd(_B, _W, TMP4);
  matAdd(TMP4, _P, _P);

  //Predict state
  ELEM(IMU_EKF_u, 0, 0) += _loop_time_s * (g_x + tmp1 * ttheta);
  ELEM(IMU_EKF_u, 1, 0) += _loop_time_s * tmp2;
  delta_u2 = _loop_time_s * (ELEM(IMU_EKF_u, 3, 0) * g_z - ELEM(IMU_EKF_u, 4, 0) * g_y - ELEM(IMU_EKF_u, 5, 0) * ELEM(IMU_EKF_u, 2, 0) - _g * stheta);
  delta_u3 = _loop_time_s * (ELEM(IMU_EKF_u, 4, 0) * g_x - ELEM(IMU_EKF_u, 2, 0) * g_z - ELEM(IMU_EKF_u, 5, 0) * ELEM(IMU_EKF_u, 3, 0) + _g * sphi * ctheta);
  ELEM(IMU_EKF_u, 2, 0) += delta_u2;
  ELEM(IMU_EKF_u, 3, 0) += delta_u3;
  ELEM(IMU_EKF_u, 4, 0) += _loop_time_s * (a_z + (_g - ELEM(IMU_EKF_u, 6, 0)) * cphi * ctheta);
  /*ELEM(IMU_EKF_u, 5,0) += 0;
   ELEM(IMU_EKF_u, 6,0) += 0;*/

  //Set Angles
  //TODO Check if they're ok
  //ahrs->e.thx = ELEM(IMU_EKF_u, 0, 0); //u(0,0) is roll according to IMU ref. frame
  //ahrs->e.thy = ELEM(IMU_EKF_u, 1, 0); //u(1,0) is pitch according to IMU ref. frame
  return;
}

//------------------------------Update with accel & gyro-------------------------------------//
void IMU_EKF_update(AHRS_State_t *ahrs, float g_x, float g_y, float g_z, float a_x, float a_y)
{
  float _Delta_m_data[2];
  Matrix Delta_m = newMatrix(2, 1, _Delta_m_data);

  //C matrix
  //_C.zeros(); //zeros or not?
  ELEM(_C, 0, 2) = -ELEM(IMU_EKF_u, 5, 0);
  ELEM(_C, 0, 3) = g_z;
  ELEM(_C, 0, 4) = -g_y;
  ELEM(_C, 0, 5) = -ELEM(IMU_EKF_u, 2, 0);
  ELEM(_C, 1, 2) = -g_z;
  ELEM(_C, 1, 3) = -ELEM(IMU_EKF_u, 5, 0);
  ELEM(_C, 1, 4) = g_x;
  ELEM(_C, 1, 5) = -ELEM(IMU_EKF_u, 3, 0);

  //Delta measures
  matZeros(Delta_m);
  ELEM(Delta_m, 0, 0) = a_x + ELEM(IMU_EKF_u, 4, 0) * g_y - ELEM(IMU_EKF_u, 3, 0) * g_z + ELEM(IMU_EKF_u, 5, 0) * ELEM(IMU_EKF_u, 2, 0);
  ELEM(Delta_m, 1, 0) = a_y + ELEM(IMU_EKF_u, 2, 0) * g_z - ELEM(IMU_EKF_u, 4, 0) * g_x + ELEM(IMU_EKF_u, 5, 0) * ELEM(IMU_EKF_u, 3, 0);

  //Gain matrix K
  //_M = QuadProd(_C, _P) + _R;
  QuadProd(_C, _P, _M);
  matAdd(_M, _R, _M);
  //_K = _P * (~_C) * (!_M);
  matMult_rhsT(_P, _C, TMP1);//TMP1 contains _P * (~_C)
  matInversed(_M, TMP2); //TMP2 contains (!_M)
  matMult(TMP1, TMP2, _K);

  //Correct state vector
  //u += _K * Delta_m
  matMult(_K, Delta_m, TMP3);
  matAdd(IMU_EKF_u, TMP3, IMU_EKF_u);

  //Updated P matrix
  //_P -= _K * _C * _P;
  //_P=(_P+(~_P))*0.5; //guarantees P to be symmetric
  matMult(_K, _C, TMP4);
  matMult(TMP4, _P, TMP5);
  matSub(_P, TMP5, _P);

  //Set Angles
  //TODO Check if they're ok
  ahrs->e.thx = ELEM(IMU_EKF_u, 0, 0); //u(0,0) is roll according to IMU ref. frame
  ahrs->e.thy = ELEM(IMU_EKF_u, 1, 0); //u(1,0) is pitch according to IMU ref. frame

  return;
}

//-------------------------Update with velocity along x,y local------------------------------//
void IMU_EKF_update_vel_xy(AHRS_State_t *ahrs, float v_x, float v_y, float dt_s)
{
  float _Delta_m_data[2];
  Matrix Delta_m = newMatrix(2, 1, _Delta_m_data);
  float _R_tmp_data[2 * 2];
  Matrix R_tmp = newMatrix(2, 2, _R_tmp_data);
  float _C_tmp_data[2 * 7];
  Matrix C_tmp = newMatrix(2, 7, _C_tmp_data);

  //R matrix
  matZeros(R_tmp);
  ELEM(R_tmp, 0, 0) = _r_vxy / dt_s;
  ELEM(R_tmp, 1, 1) = ELEM(R_tmp, 0, 0);

  //C matrix
  matZeros(C_tmp);
  ELEM(C_tmp, 0, 2) = 1.f;
  ELEM(C_tmp, 1, 3) = 1.f;

  //Delta measures
  matZeros(Delta_m);
  ELEM(Delta_m, 0, 0) = v_x - ELEM(IMU_EKF_u, 2, 0);
  ELEM(Delta_m, 1, 0) = v_y - ELEM(IMU_EKF_u, 3, 0);

  //Gain matrix K
  //_M = QuadProd(C_tmp,_P) + R_tmp;
  QuadProd(C_tmp, _P, _M);
  matAdd(_M, R_tmp, _M);
  //_K = _P * (~C_tmp) * (!_M);
  matMult_rhsT(_P, C_tmp, TMP1); //TMP1 contains _P * (~_C)
  matInversed(_M, TMP2); //TMP2 contains (!_M)
  matMult(TMP1, TMP2, _K);

  //Correct state vector
  //u += _K * Delta_m;
  matMult(_K, Delta_m, TMP3);
  matAdd(IMU_EKF_u, TMP3, IMU_EKF_u);

  //Updated P matrix
  //_P -= _K * C_tmp * _P;
  //_P=(_P+(~_P))*0.5; //guarantees P to be symmetric
  matMult(_K, C_tmp, TMP4);
  matMult(TMP4, _P, TMP5);
  matSub(_P, TMP5, _P);

  //Set Angles
  //TODO Check if they're ok
  ahrs->e.thx = ELEM(IMU_EKF_u, 0, 0); //u(0,0) is roll according to IMU ref. frame
  ahrs->e.thy = ELEM(IMU_EKF_u, 1, 0); //u(1,0) is pitch according to IMU ref. frame

  return;
}

//--------------------------Update with velocity along z local--------------------------------//
void IMU_EKF_update_vel_z(AHRS_State_t *ahrs, float v_z, float dt_s)
{
  float _C_tmp_data[1 * 7];
  Matrix C_tmp = newMatrix(1, 7, _C_tmp_data);
  float _K_tmp_data[7 * 1];
  Matrix K = newMatrix(7, 1, _K_tmp_data);

  //C matrix
  matZeros(C_tmp);
  ELEM(C_tmp, 0, 4) = 1.f;

  //Delta measures
  float Delta_m = v_z - ELEM(IMU_EKF_u, 4, 0);

  //Gain matrix K (change K with _K in correction and update of P)
  /*_M = QuadProd(C_tmp,_P) + (_r_vz / dt_s);
   _K = _P * (~C_tmp) * (!_M);*/

  //Faster Gain matrix K
  float inv_m = 1.f / (ELEM(_P, 4, 4) + (_r_vz / dt_s));
  ELEM(K, 0, 0) = ELEM(_P, 0, 4) * inv_m;
  ELEM(K, 1, 0) = ELEM(_P, 1, 4) * inv_m;
  ELEM(K, 2, 0) = ELEM(_P, 2, 4) * inv_m;
  ELEM(K, 3, 0) = ELEM(_P, 3, 4) * inv_m;
  ELEM(K, 4, 0) = ELEM(_P, 4, 4) * inv_m;
  ELEM(K, 5, 0) = ELEM(_P, 5, 4) * inv_m;
  ELEM(K, 6, 0) = ELEM(_P, 6, 4) * inv_m;

  //Correct state vector
  //u += K * Delta_m;
  matMultScalar(K, Delta_m, TMP3);
  matAdd(IMU_EKF_u, TMP3, IMU_EKF_u);

  //Updated P matrix
  //_P -= K * C_tmp * _P;
  //_P=(_P+(~_P))*0.5; //guarantees P to be symmetric
  matMult(K, C_tmp, TMP4);
  matMult(TMP4, _P, TMP5);
  matSub(_P, TMP5, _P);

  //Set Angles
  //TODO Check if they're ok
  ahrs->e.thx = ELEM(IMU_EKF_u, 0, 0); //u(0,0) is roll according to IMU ref. frame
  ahrs->e.thy = ELEM(IMU_EKF_u, 1, 0); //u(1,0) is pitch according to IMU ref. frame

  return;
}

//--------------------------Update with velocity along d global--------------------------------//
void IMU_EKF_update_vel_d(AHRS_State_t *ahrs, float v_d, float dt_s)
{
  float _C_tmp_data[1 * 7];
  Matrix C_tmp = newMatrix(1, 7, _C_tmp_data);
  float _K_tmp_data[7 * 1];
  Matrix K = newMatrix(7, 1, _K_tmp_data);
  float _M_tmp_data[1 * 1];
  Matrix M = newMatrix(1, 1, _M_tmp_data);

  //Trig functions
  float sphi = sinf(ELEM(IMU_EKF_u, 0, 0));
  float cphi = cosf(ELEM(IMU_EKF_u, 0, 0));
  float stheta = sinf(ELEM(IMU_EKF_u, 1, 0));
  float ctheta = cosf(ELEM(IMU_EKF_u, 1, 0));

  //C matrix
  ELEM(C_tmp, 0, 0) = ELEM(IMU_EKF_u, 3,0) * cphi * ctheta - ELEM(IMU_EKF_u, 4,0) * ctheta * sphi;
  ELEM(C_tmp, 0, 1) = -ELEM(IMU_EKF_u, 2, 0) * ctheta - ELEM(IMU_EKF_u, 4,0) * cphi * stheta - ELEM(IMU_EKF_u, 3,0) * sphi * stheta;
  ELEM(C_tmp, 0, 2) = -stheta;
  ELEM(C_tmp, 0, 3) = ctheta * sphi;
  ELEM(C_tmp, 0, 4) = cphi * ctheta;
  ELEM(C_tmp, 0, 5) = 0;
  ELEM(C_tmp, 0, 6) = 0;

  //Delta measures
  float Delta_m = v_d - (ELEM(IMU_EKF_u, 4,0) * cphi * ctheta - ELEM(IMU_EKF_u, 2,0) * stheta + ELEM(IMU_EKF_u, 3,0) * ctheta * sphi);

  //Gain matrix K
  //_M = QuadProd(C_tmp,_P) + (_r_vd / dt_s);
  QuadProd(C_tmp, _P, M);
  ELEM(M, 0 , 0) += (_r_vd / dt_s);
  //K = _P * (~C_tmp) * (!_M);
  matMult_rhsT(_P, C_tmp, K);
  matMultScalar(K, 1.0f / ELEM(M, 0, 0), K);

  //Correct state vector
  //u += K * Delta_m;
  matMultScalar(K, Delta_m, TMP3);
  matAdd(IMU_EKF_u, TMP3, IMU_EKF_u);

  //Updated P matrix
  //_P -= _K * C_tmp * _P;
  //_P=(_P+(~_P))*0.5; //guarantees P to be symmetric
  matMult(K, C_tmp, TMP4);
  matMult(TMP4, _P, TMP5);
  matSub(_P, TMP5, _P);

  //Set Angles
  //TODO Check if they're ok
  ahrs->e.thx = ELEM(IMU_EKF_u, 0, 0); //u(0,0) is roll according to IMU ref. frame
  ahrs->e.thy = ELEM(IMU_EKF_u, 1, 0); //u(1,0) is pitch according to IMU ref. frame

  return;
}

//-----------------------------------Starting values-----------------------------------------//
void IMU_EKF_set_starting_values(AHRS_State_t *ahrs, float Phi_0, float Theta_0, float c_damp_0)
{
  ELEM(IMU_EKF_u, 0,0) = Phi_0;
  ELEM(IMU_EKF_u, 1,0) = Theta_0;
  ELEM(IMU_EKF_u, 5,0) = c_damp_0;

  //Set Angles
  //TODO Check if they're ok
  ahrs->e.thx = ELEM(IMU_EKF_u, 0, 0); //u(0,0) is roll according to IMU ref. frame
  ahrs->e.thy = ELEM(IMU_EKF_u, 1, 0); //u(1,0) is pitch according to IMU ref. frame

  return;
}

//-------------------------------------Input noises------------------------------------------//
void IMU_EKF_set_input_noises(float g_xy, float g_z, float a_z, float c_damp, float b_az)
{
  ELEM(_W, 0, 0) = g_xy * _loop_time_s;
  ELEM(_W, 1, 1) = g_xy * _loop_time_s;
  ELEM(_W, 2, 2) = g_z * _loop_time_s;
  ELEM(_W, 3, 3) = a_z * _loop_time_s;
  ELEM(_W, 4, 4) = c_damp * _loop_time_s;
  ELEM(_W, 5, 5) = b_az * _loop_time_s;

  return;
}

//----------------------------------Output accel noises--------------------------------------/
void IMU_EKF_set_acc_noise(float a)
{
  float inv_loop_time_s = 1.0f / _loop_time_s;
  ELEM(_R, 0, 0) = a * inv_loop_time_s;
  ELEM(_R, 1, 1) = a * inv_loop_time_s;

  return;
}

//-------------------------Output velocity noises along x,y local-----------------------------/
void IMU_EKF_set_vel_xy_noise(float v_xy)
{
  _r_vxy = v_xy;

  return;
}

//--------------------------Output velocity noises along z local------------------------------/
void IMU_EKF_set_vel_z_noise(float v_z)
{
  _r_vz = v_z;

  return;
}

//-------------------------Output velocity noises along d global------------------------------/
void IMU_EKF_set_vel_d_noise(float v_d)
{
  _r_vd = v_d;

  return;
}
