//
//  IMU_EKF_AV.h
//
//
//  Created by Andrea Vivani on 25/06/22.
//  Copyright (c) 2022 Andrea Vivani. All rights reserved.
//
// AXES DIRECTIONS (SENSORS AND ESTIMATION): X POINTING FORWARD (ROLL, PHI), Y POINTING RIGHT (PITCH, THETA), Z POINTING DOWN (YAW, PSI)
// IF SENSORS HAVE DIFFERENT ORIENTATION, ROTATE READINGS PRIOR THAN INPUT THEM TO THE ESTIMATORS

#ifndef IMU_EKF_AV_H_
#define IMU_EKF_AV_H_
#include "Matrix.h"
#include "math.h"
#include "ahrs.h"

//======================================Parameters=============================================//
//----------------------------------------Noises-----------------------------------------------//
#define IMU_EKF_AV_GXY_NOISE 	1e-3f    //gyro x,y noise
#define IMU_EKF_AV_GZ_NOISE 	1e-3f    //gyro z noise
#define IMU_EKF_AV_AXY_NOISE    1e-2f	//accel x,y noise
#define IMU_EKF_AV_AZ_NOISE  	1e-1f	//accel z noise
#define IMU_EKF_AV_C_DAMP_NOISE	5e-4f	//damping coefficient noise
#define IMU_EKF_AV_B_AZ_NOISE	1e-2f	//bias acc z noise
#define IMU_EKF_AV_VXY_NOISE    1e-3f   //velocity x,y local noise
#define IMU_EKF_AV_VZ_NOISE     1e-3f   //velocity z local noise
#define IMU_EKF_AV_VD_NOISE     1e-3f   //velocity d global noise (pointing down)

//------------------------------------Initial Values-------------------------------------------//
#define IMU_EKF_AV_PHI0         0.f		//initial value of the Roll angle, in rad
#define IMU_EKF_AV_THETA0       0.f		//initial value of the Pitch angle, in rad
#define IMU_EKF_AV_C_DAMP0      0.07413f	//initial value of the damping coefficient, inversely proportional to the mass of the quadrotor, in N*s/m

//Roll and Pitch to be part of ahrs.e.thx/thy (Euler angles in RPY order: Phi=roll, Teta=pitch, Psi=yaw), in rad
//Translational velocities are ahrs.Xd/Yd/Zd

void IMU_EKF_init(AHRS_State_t *ahrs, float g_val, float loop_time_s);
void IMU_EKF_prediction(float g_x, float g_y, float g_z, float a_z);
void IMU_EKF_update(AHRS_State_t *ahrs, float g_x, float g_y, float g_z, float a_x, float a_y);
void IMU_EKF_update_vel_xy(AHRS_State_t *ahrs, float v_x, float v_y, float dt_s);
void IMU_EKF_update_vel_z(AHRS_State_t *ahrs, float v_z, float dt_s);
void IMU_EKF_update_vel_d(AHRS_State_t *ahrs, float v_d, float dt_s);
void IMU_EKF_set_starting_values(AHRS_State_t *ahrs, float Phi_0, float Theta_0, float c_damp_0);
void IMU_EKF_set_input_noises(float g_xy, float g_z, float a_z, float c_damp, float b_az);
void IMU_EKF_set_acc_noise(float a);
void IMU_EKF_set_vel_xy_noise(float v_xy);
void IMU_EKF_set_vel_z_noise(float v_z);
void IMU_EKF_set_vel_d_noise(float v_d);

#endif

