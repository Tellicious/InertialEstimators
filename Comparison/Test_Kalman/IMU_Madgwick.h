//  IMU_Madgwick.h
//
//
//  Created by Sebastian Madgwick
//
//
#ifndef IMU_Madgwick_H_
#define IMU_Madgwick_H_
#include "math.h"

//======================================Parameters=============================================//
#define IMU_MADGWICK_GYRO_ERROR 0.0872664626f //gyroscope supposed measurement error (rad/s)
#define IMU_MADGWICK_GYRO_DRIFT 0.0034906585f //gyroscope supposed bias drift (rad/s/s)

class IMU_Madgwick{
	public:
		float Roll, Pitch; //estimated orientation (Euler angles in RPY order: Phi=Roll, Teta=Pitch), in rad
		float q1, q2, q3, q4; //estimated quaternion
		IMU_Madgwick(float gyro_error=IMU_MADGWICK_GYRO_ERROR);
		void compute(float g_x, float g_y, float g_z, float a_x, float a_y, float a_z, float loop_time_s);	//compute estimation
	private:
		float _beta;	//parameter depending on supposed measurement error
};
#endif