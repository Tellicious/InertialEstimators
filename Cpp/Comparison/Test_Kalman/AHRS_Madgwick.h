//  AHRS_Madgwick.h
//
//
//  Created by Sebastian Madgwick
//
//
#ifndef AHRS_Madgwick_H_
#define AHRS_Madgwick_H_

#include "math.h"

//======================================Parameters=============================================//
#define AHRS_MADGWICK_GYRO_ERROR 0.0872664626f //gyroscope supposed measurement error (rad/s)
#define AHRS_MADGWICK_GYRO_DRIFT 0.0034906585f //gyroscope supposed bias drift (rad/s/s)

class AHRS_Madgwick{
	public:
		float Roll, Pitch, Yaw; //estimated orientation (Euler angles in RPY order: Phi=Roll, Teta=Pitch, Psi=Yaw), in rad
		float q1, q2, q3, q4; //estimated quaternion
		AHRS_Madgwick(float gyro_error=AHRS_MADGWICK_GYRO_ERROR, float gyro_drift=AHRS_MADGWICK_GYRO_DRIFT);
		void compute(float g_x, float g_y, float g_z, float a_x, float a_y, float a_z, float m_x, float m_y, float m_z, float loop_time_s);	//compute estimation
	private:
		float _beta;	//parameter depending on supposed measurement error
		float _zeta;	//parameter depending on supposed drift
		float _b_x, _b_z; 	//magnetic field components
		float _w_bx, _w_by, _w_bz; //estimated gyroscope biases
};
#endif