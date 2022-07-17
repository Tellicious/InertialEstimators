#ifndef _ALTITUDE_KF_H
#define _ALTITUDE_KF_H
#include <stm32f10x.h>

void alt_KF_init (uint32_t loop_time_ms);
void alt_KF_prediction (void);
void alt_KF_update_baro (void);
void alt_KF_update_accel (void);
void alt_KF_update_INS(void);

#endif //_ALTITUDE_KF_H
