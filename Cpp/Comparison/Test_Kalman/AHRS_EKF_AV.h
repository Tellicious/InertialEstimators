//  AHRS_EKF_AV.h
//
//
//  Created by Andrea Vivani on 21/4/15.
//  Copyright (c) 2015 Andrea Vivani. All rights reserved.
//
#ifndef AHRS_EKF_AV_H_
#define AHRS_EKF_AV_H_
#include "MatLib.h"
#include "math.h"

//======================================Parameters=============================================//
//----------------------------------------Noises-----------------------------------------------//
#define AHRS_EKF_AV_GXY_NOISE       1e-3f	//gyro x,y noise
#define AHRS_EKF_AV_GZ_NOISE        1e-3f	//gyro z noise
#define AHRS_EKF_AV_A_NOISE         1e-2f	//accel x and y noise
#define AHRS_EKF_AV_AZ_NOISE        1e-1f	//accel z noise
#define AHRS_EKF_AV_M_NOISE         1e-2f	//mag noise
#define AHRS_EKF_AV_C_DAMP_NOISE	5e-4f	//damping coefficient noise
#define AHRS_EKF_AV_B_AZ_NOISE      1e-3f	//bias acc z noise
#define AHRS_EKF_AV_INCL_NOISE      1e-3f	//bias magnetic field inclination noise
#define AHRS_EKF_AV_B_G_NOISE       1e-4f	//bias gyro x noise
#define AHRS_EKF_AV_VXY_NOISE       1e-3f   //velocity x,y local noise
#define AHRS_EKF_AV_VZ_NOISE        1e-3f   //velocity z local noise
#define AHRS_EKF_AV_VNE_NOISE       1e-3f   //velocity n,e global noise
#define AHRS_EKF_AV_VD_NOISE        1e-3f   //velocity d global noise (pointing down)

//------------------------------------Initial Values-------------------------------------------//
#define AHRS_EKF_AV_PHI0            0.f		//initial value of the Roll angle, in rad
#define AHRS_EKF_AV_THETA0          0.f		//initial value of the Pitch angle, in rad
#define AHRS_EKF_AV_PSI0            0.f		//initial value of the Yaw angle, in rad
#define AHRS_EKF_AV_C_DAMP0         0.07413f	//initial value of the damping coefficient, inversely proportional to the mass of the quadrotor, in N*s/m
#define AHRS_EKF_AV_INCL0           -1.0734f	//initial value of the inclination of the magnetic field, in rad (Northern Emisphere positive (pointing down), Southern Emisphere negative (pointing up))

class AHRS_EKF_AV{
    
	public:
		float Roll, Pitch, Yaw; //estimated orientation (Euler angles in RPY order: Phi=roll, Teta=pitch, Psi=yaw), in rad
		float Xd, Yd, Zd; //estimated translational velocities in local frame of reference, in m/s
		MatrixXf u;	// state vector (Roll=Phi, Pitch=Theta, Yaw=Psi, Xd, Yd, Zd, c_damp, b_az, incl, b_gx, b_gy, b_gz) angles in rad, angular velocities in rad/s, velocities in m/s, c_damp in N*s/m
		AHRS_EKF_AV(float g_val, float loop_time_pred_update_acc_s, float loop_time_mag_s);
        AHRS_EKF_AV(float g_val, float loop_time_pred_update_acc_s);
		void prediction(float g_x, float g_y, float g_z, float a_z);
		void update_accel(float g_x, float g_y, float g_z, float a_x, float a_y);
		void update_mag(float m_x, float m_y, float m_z);
        void update_vel_xy(float v_x, float v_y, float dt_s);
        void update_vel_z(float v_z, float dt_s);
        void update_vel_ne(float v_n, float v_e, float dt_s);
        void update_vel_d(float v_d, float dt_s);
		void set_starting_values(float Phi_0=AHRS_EKF_AV_PHI0, float Theta_0=AHRS_EKF_AV_THETA0, float Psi_0=AHRS_EKF_AV_PSI0, float c_damp_0=AHRS_EKF_AV_C_DAMP0,
			float incl_0=AHRS_EKF_AV_INCL0);
		void set_input_noises(float g_xy=AHRS_EKF_AV_GXY_NOISE, float g_z=AHRS_EKF_AV_GZ_NOISE, float a_z=AHRS_EKF_AV_AZ_NOISE, float c_damp=AHRS_EKF_AV_C_DAMP_NOISE, float b_az=AHRS_EKF_AV_B_AZ_NOISE, float incl=AHRS_EKF_AV_INCL_NOISE, float b_g=AHRS_EKF_AV_B_G_NOISE);
		void set_acc_noises(float a=AHRS_EKF_AV_A_NOISE);
		void set_mag_noises(float m=AHRS_EKF_AV_M_NOISE);
        void set_vel_xy_noise(float v_xy=AHRS_EKF_AV_VXY_NOISE){_r_vxy=v_xy;return;};
        void set_vel_z_noise(float v_z=AHRS_EKF_AV_VZ_NOISE){_r_vz=v_z;return;};
        void set_vel_ne_noise(float v_ne=AHRS_EKF_AV_VNE_NOISE){_r_vne=v_ne;return;};
        void set_vel_d_noise(float v_d=AHRS_EKF_AV_VD_NOISE){_r_vd=v_d;return;};
    
	private:
		MatrixXf _A;        //state matrix
		MatrixXf _B;        //input matrix
		MatrixXf _C_acc;	//acceleration output matrix
		MatrixXf _C_mag;	//magnetic field output matrix
		MatrixXf _P;        //errors expected value matrix
		MatrixXf _W;        //gyro and acc_z noises covariance matrix (g_x, g_y, g_z, a_z, c_damp, b_az, incl, b_gx, b_gy, b_gz)
		MatrixXf _R_acc; 	//acc_x and acc_y noises covariance matrix
		MatrixXf _R_mag; 	//magnetometer noises covariance matrix
		MatrixXf _M;		//temporary matrix
		MatrixXf _K;		//gain matrix
		float _loop_time_pred_update_acc_s;	//prediction and update with accelerometer loop time, in s
		float _loop_time_mag_s;             //update with magnetometer loop time, in s
		float _g;           //value of gravitational acceleration, in m/s^2
        float _r_vxy, _r_vz, _r_vne, _r_vd; //velocities noise covariances
};
#endif