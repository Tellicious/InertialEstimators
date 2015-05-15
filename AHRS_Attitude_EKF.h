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
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 ****************************************************************************/

/*
 * @file AHRS_Attitude_EKF_H_
 *
 * Extended Kalman Filter for Attitude Estimation.
 *
 * @author Tobias Naegeli <naegelit@student.ethz.ch>
 * @author Lorenz Meier <lm@inf.ethz.ch>
 * @author Thomas Gubler <thomasgubler@gmail.com>
 */

#ifndef AHRS_ATTITUDE_EKF_H_
#define AHRS_ATTITUDE_EKF_H_
#include "MatLib.h"
#include "math.h"

//#define AHRS_ATTITUDE_EKF_USE_COMPLETE_INERTIA_MATRIX
//#define AHRS_ATTITUDE_EKF_USE_DIAGONAL_INERTIA_MATRIX

//======================================Parameters=============================================//
//----------------------------------------Noises-----------------------------------------------//
#define AHRS_ATTITUDE_EKF_R_S_NOISE 	2e-2f	//rotational speed process noise
#define AHRS_ATTITUDE_EKF_R_A_NOISE     16.f	//rotational acceleration process noise
#define AHRS_ATTITUDE_EKF_A_NOISE       1.8f	//acceleration process noise
#define AHRS_ATTITUDE_EKF_M_NOISE       1.f	//magnetic field process noise
#define AHRS_ATTITUDE_EKF_GYRO_NOISE	4e-6f	//gyroscope noise
#define AHRS_ATTITUDE_EKF_ACCEL_NOISE	50.f	//accelerometer noise
#define AHRS_ATTITUDE_EKF_MAG_NOISE     0.5f	//magnetometer noise

//------------------------------------Initial Values-------------------------------------------//
#define AHRS_ATTITUDE_EKF_INCL0         -1.0734f	//initial value of the inclination of the magnetic field, in rad (Northern Emisphere positive (pointing down), Southern Emisphere negative (pointing up))


class AHRS_Attitude_EKF {
public:
    float Roll, Pitch, Yaw; //estimated orientation (Euler angles in RPY order: Phi=roll, Teta=pitch, Psi=yaw), in rad
    MatrixXf u; //estimated state (3 angular velocities, in rad/s, 3 angular accelerations, in rad/s^2, acceleration components in body frame, in m/s^2, magnetic field components)
    AHRS_Attitude_EKF(float g_val, float loop_time_s);
    AHRS_Attitude_EKF(float g_val, float loop_time_s, MatrixXf &inertia_m);
    AHRS_Attitude_EKF(float g_val, float loop_time_pred_s, float loop_time_update_gyro_s, float loop_time_update_acc_s, float loop_time_update_mag_s, MatrixXf &inertia_m);
    void prediction();
    void update_gyro(float g_x, float g_y, float g_z);
    void update_accel(float a_x, float a_y, float a_z);
    void update_mag(float m_x, float m_y, float m_z);
    void euler_angles(); //compute Euler angles from state estimation
    void set_inclination(float incl_angle=AHRS_ATTITUDE_EKF_INCL0); //set magnetic field inclination angle
    void set_process_noise(float rot_sp=AHRS_ATTITUDE_EKF_R_S_NOISE, float rot_acc=AHRS_ATTITUDE_EKF_R_A_NOISE, float acc=AHRS_ATTITUDE_EKF_A_NOISE, float mag=AHRS_ATTITUDE_EKF_M_NOISE);
    void set_gyro_noise(float g=AHRS_ATTITUDE_EKF_GYRO_NOISE);
    void set_accel_noise(float a=AHRS_ATTITUDE_EKF_ACCEL_NOISE);
    void set_mag_noise(float m=AHRS_ATTITUDE_EKF_MAG_NOISE);
    
private:
    void initialize(float g_val, float loop_time_pred_s, float loop_time_update_gyro_s, float loop_time_update_acc_s, float loop_time_update_mag_s, MatrixXf inertia_m); //initialization function
    MatrixXf _A; //state matrix
    MatrixXf _C; //variable output matrix
    MatrixXf _J; //inertia matrix
    MatrixXf _Ji; //inverse of inertia matrix
    MatrixXf _Q;  //process noise covariance matrix
    MatrixXf _M; //temp matrix
    MatrixXf _K; //gain matrix
    MatrixXf _R_gyro; //gyro noise covariance matrix
    MatrixXf _R_accel; //accel noise covariance matrix
    MatrixXf _R_mag;    //mag noise covariance matrix
    MatrixXf _P;    //expected error values matrix
    
    float _loop_time_pred_s; //prediction loop time, in s
    float _loop_time_update_gyro_s; //update with gyro loop time, in s
    float _loop_time_update_acc_s;  //update with acc loop time, in s
    float _loop_time_update_mag_s;  //update with mag loop time, in s
    float _g_val; //value of gravitrational acceleration, in m/s^s
};

#endif