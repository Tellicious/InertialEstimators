//  AHRS_Madgwick.cpp
//
//
//  Created by Sebastian Madgwick
//
//
#include "AHRS_Madgwick.h"

//=====================================Constructor===========================================//
AHRS_Madgwick::AHRS_Madgwick(float gyro_error, float gyro_drift){
	// Initial quaternion
	q1=1;
	q2=0;
	q3=0;
	q4=0;
	// Initial magnetic field estimated direction
	_b_x=1;
	_b_z=0;
	// Initial gyroscope estimated biases
	_w_bx=0;
	_w_by=0;
	_w_bz=0;
	// filter parameters
	_beta=sqrtf(3.0f*0.25f)*gyro_error;
	_zeta=sqrtf(3.0f*0.35f)*gyro_drift;
}

//====================================Public Members==========================================//
//--------------------------Compute---------------------------//
void AHRS_Madgwick::compute(float g_x, float g_y, float g_z, float a_x, float a_y, float a_z, float m_x, float m_y, float m_z, float loop_time_s){
	// local system variables
	float inv_norm;	//inverse of vector norm
	float SEqDot_omega_1, SEqDot_omega_2, SEqDot_omega_3, SEqDot_omega_4; //quaternion rate from gyroscopes elements
	float f_1, f_2, f_3, f_4, f_5, f_6;	//objective function elements
	float J_11or24, J_12or23, J_13or22, J_14or21, J_32, J_33,	//objective function Jacobian elements
	J_41, J_42, J_43, J_44, J_51, J_52, J_53, J_54, J_61, J_62, J_63, J_64;
	float SEqHatDot_1, SEqHatDot_2, SEqHatDot_3, SEqHatDot_4;//estimated direction of the gyroscope error
	float w_err_x, w_err_y, w_err_z; //estimated direction of the gyroscope error (angular)
	float h_x, h_y, h_z;	//computed flux in the earth frame

	// axulirary variables to avoid reapeated calculations
	float halfSEq_1 = 0.5f * q1;
	float halfSEq_2 = 0.5f * q2;
	float halfSEq_3 = 0.5f * q3;
	float halfSEq_4 = 0.5f * q4; 
	float twoSEq_1 = 2.0f * q1;
	float twoSEq_2 = 2.0f * q2; 
	float twoSEq_3 = 2.0f * q3; 
	float twoSEq_4 = 2.0f * q4; 
	float twob_x = 2.0f * _b_x; 
	float twob_z = 2.0f * _b_z;
	float twob_xSEq_1 = 2.0f * _b_x * q1; 
	float twob_xSEq_2 = 2.0f * _b_x * q2; 
	float twob_xSEq_3 = 2.0f * _b_x * q3; 
	float twob_xSEq_4 = 2.0f * _b_x * q4; 
	float twob_zSEq_1 = 2.0f * _b_z * q1; 
	float twob_zSEq_2 = 2.0f * _b_z * q2; 
	float twob_zSEq_3 = 2.0f * _b_z * q3; 
	float twob_zSEq_4 = 2.0f * _b_z * q4; 
	float SEq_1SEq_2;
	float SEq_1SEq_3 = q1 * q3; 
	float SEq_1SEq_4;
	float SEq_2SEq_3;
	float SEq_2SEq_4 = q2 * q4; 
	float SEq_3SEq_4;
    float SEq_1SEq_1;
    float SEq_2SEq_2;
    float SEq_3SEq_3;
    float SEq_4SEq_4;
	float twom_x = 2.0f * m_x;
	float twom_y = 2.0f * m_y;
	float twom_z = 2.0f * m_z;

	// normalize the accelerometer measurement 
	inv_norm = -1.f/sqrtf(a_x * a_x + a_y * a_y + a_z * a_z);
    if (isnan(inv_norm)||isinf(inv_norm)){
        inv_norm=1.f;
    }
	a_x *= inv_norm;
	a_y *= inv_norm;
	a_z *= inv_norm;
    
	// normalize the magnetometer measurement
	inv_norm = 1.f/sqrtf(m_x * m_x + m_y * m_y + m_z * m_z);
    if (isnan(inv_norm)||isinf(inv_norm)){
        inv_norm=1.f;
    }
	m_x *= inv_norm;
	m_y *= inv_norm;
	m_z *= inv_norm;

	// compute the objective function and Jacobian
	f_1 = twoSEq_2 * q4 - twoSEq_1 * q3 - a_x;
	f_2 = twoSEq_1 * q2 + twoSEq_3 * q4 - a_y;
	f_3 = 1.0f - twoSEq_2 * q2 - twoSEq_3 * q3 - a_z;
	f_4 = twob_x * (0.5f - q3 * q3 - q4 * q4) + twob_z * (SEq_2SEq_4 - SEq_1SEq_3) - m_x; 
	f_5 = twob_x * (q2 * q3 - q1 * q4) + twob_z * (q1 * q2 + q3 * q4) - m_y;
	f_6 = twob_x * (SEq_1SEq_3 + SEq_2SEq_4) + twob_z * (0.5f - q2 * q2 - q3 * q3) - m_z;
	J_11or24 = twoSEq_3;
	J_12or23 = 2.0f * q4;
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
	inv_norm = 1.f/sqrtf(SEqHatDot_1 * SEqHatDot_1 + SEqHatDot_2 * SEqHatDot_2 + SEqHatDot_3 * SEqHatDot_3 + SEqHatDot_4 * SEqHatDot_4);
    if (isnan(inv_norm)||isinf(inv_norm)){
        inv_norm=1.f;
    }
	SEqHatDot_1 = SEqHatDot_1 * inv_norm;
	SEqHatDot_2 = SEqHatDot_2 * inv_norm;
	SEqHatDot_3 = SEqHatDot_3 * inv_norm;
	SEqHatDot_4 = SEqHatDot_4 * inv_norm;

	// compute angular estimated direction of the gyroscope error
	w_err_x = twoSEq_1 * SEqHatDot_2 - twoSEq_2 * SEqHatDot_1 - twoSEq_3 * SEqHatDot_4 + twoSEq_4 * SEqHatDot_3; 
	w_err_y = twoSEq_1 * SEqHatDot_3 + twoSEq_2 * SEqHatDot_4 - twoSEq_3 * SEqHatDot_1 - twoSEq_4 * SEqHatDot_2; 
	w_err_z = twoSEq_1 * SEqHatDot_4 - twoSEq_2 * SEqHatDot_3 + twoSEq_3 * SEqHatDot_2 - twoSEq_4 * SEqHatDot_1;

	// compute and remove the gyroscope biases 
	_w_bx += w_err_x * loop_time_s * _zeta;
	_w_by += w_err_y * loop_time_s * _zeta;
	_w_bz += w_err_z * loop_time_s * _zeta;
	g_x -= _w_bx; 
	g_y -= _w_by; 
	g_z -= _w_bz;

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
	inv_norm = 1.f/sqrtf((q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4));
    if (isnan(inv_norm)||isinf(inv_norm)){
        inv_norm=1.f;
    }
	q1 *= inv_norm;
	q2 *= inv_norm;
	q3 *= inv_norm;
	q4 *= inv_norm;

	// compute flux in the earth frame 
	SEq_1SEq_2 = q1 * q2; 
	SEq_1SEq_3 = q1 * q3; 
	SEq_1SEq_4 = q1 * q4; 
	SEq_3SEq_4 = q3 * q4; 
	SEq_2SEq_3 = q2 * q3; 
	SEq_2SEq_4 = q2 * q4;
    SEq_1SEq_1 = q1 * q1;
    SEq_2SEq_2 = q2 * q2;
    SEq_3SEq_3 = q3 * q3;
    SEq_4SEq_4 = q4 * q4;
	h_x = twom_x * (0.5f - SEq_3SEq_3 - SEq_4SEq_4) + twom_y * (SEq_2SEq_3 - SEq_1SEq_4) + twom_z * (SEq_2SEq_4 + SEq_1SEq_3);
	h_y = twom_x * (SEq_2SEq_3 + SEq_1SEq_4) + twom_y * (0.5f - SEq_2SEq_2 - SEq_4SEq_4) + twom_z * (SEq_3SEq_4 - SEq_1SEq_2);
	h_z = twom_x * (SEq_2SEq_4 - SEq_1SEq_3) + twom_y * (SEq_3SEq_4 + SEq_1SEq_2) + twom_z * (0.5f - SEq_2SEq_2 - SEq_3SEq_3);

	// normalize the flux vector to have only components in the x and z 
	_b_x = sqrtf((h_x * h_x) + (h_y * h_y));
	_b_z = h_z;

	//convert from quaternion to Euler angles
    Roll = atan2f((2.f * (SEq_3SEq_4 + SEq_1SEq_2)), (SEq_1SEq_1 - SEq_2SEq_2 - SEq_3SEq_3 + SEq_4SEq_4));
    Pitch = -asinf((2.f * (SEq_2SEq_4 - SEq_1SEq_3)));
    Yaw = atan2f((2.f * (SEq_2SEq_3 + SEq_1SEq_4)), (SEq_1SEq_1 + SEq_2SEq_2 - SEq_3SEq_3 - SEq_4SEq_4));
	return;
}