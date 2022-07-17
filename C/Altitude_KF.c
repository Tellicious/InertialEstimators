#include "Altitude_KF.h"
#include "GlobalVariables.h"
#include "ConfigParams.h"
#include <math.h>

typedef struct{
	float h, a, v;
	float m_hh, m_ha, m_ah, m_aa, m_vh, m_va; //gain matrix
	float T, T2_2; //step time
} alt_KF_t;

alt_KF_t alt_KF;

/* MATLAB CODE
T = 0.01;
A=[1 T T^2/2; 0 1 T; 0 0 1];
B = eye(3);
C= [1 0 0; 0 0 1];
Plant = ss(A,B,C,0,T,'inputname',{'u' 'v' 'w'},'outputname',{'h' 'a'});
Q = [10 0 0; 0 100 0; 0 0 100];
R = [10000 0; 0 1000000];
[kalmf,L,P,M] = kalman(Plant,Q,R);
M
*/
void alt_KF_init (uint32_t loop_time_ms){
	alt_KF.T = loop_time_ms * 0.001;
	alt_KF.T2_2 = (alt_KF.T * alt_KF.T) / 2;
	alt_KF.m_hh = 0.057110675473981;
	alt_KF.m_vh = 0.120442911328370;
	alt_KF.m_ah = 0.046044355080117;
	alt_KF.m_ha = 0.000460443550801;
	alt_KF.m_va = 0.003177022298499;
	alt_KF.m_aa = 0.008743105577194;
	return;
}

void alt_KF_prediction (void){
	alt_KF.h += alt_KF.T2_2 * alt_KF.a + alt_KF.T * alt_KF.v;
	alt_KF.v += alt_KF.T * alt_KF.a;
	return;
}












void alt_KF_update_baro (void){
	float alt_delta = 44330.76923 * (1 - powf((INSData.pressure * NAVData.inv_referencePressure), 0.190266436)) - alt_KF.h;
	alt_KF.h += alt_KF.m_hh * alt_delta;
	alt_KF.a += alt_KF.m_ah * alt_delta;
	alt_KF.v += alt_KF.m_vh * alt_delta;
	alt_KF_update_INS();
	return;
}

void alt_KF_update_accel (void){
	float accel_delta = INS_G_VAL - INSData.cal_accel[2] - alt_KF.a; //accel sign to be checked
	alt_KF.h += alt_KF.m_ha * accel_delta;
	alt_KF.a += alt_KF.m_aa * accel_delta;
	alt_KF.v += alt_KF.m_va * accel_delta;
	alt_KF_update_INS();
	return;
}

void alt_KF_update_INS(void){
	NAVData.altitude = alt_KF.h;
	NAVData.rateOfClimb = alt_KF.v;
	return;
}
