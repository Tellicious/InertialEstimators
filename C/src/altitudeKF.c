/* BEGIN Header */
/**
 ******************************************************************************
 * \file            altitudeKF.c
 * \author          Andrea Vivani
 * \brief           Kalman filter for altitude estimation
 ******************************************************************************
 * \copyright
 *
 * Copyright 2023 Andrea Vivani
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the “Software”), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 *
 ******************************************************************************
 */
/* END Header */

/* Includes ------------------------------------------------------------------*/

#include "altitudeKF.h"
#include <math.h>
#include "basicMath.h"
#include "string.h"

#if defined(configALTITUDE_KF_ACC_HP_FILTER) || (configUSE_ALT_TOF != configTOF_DISABLE)
#include "IIRfilters.h"
#endif

/* MATLAB Code ---------------------------------------------------------------*/

/* MATLAB CODE for KF
 * T=1/200;
 * A=[1 T T^2/2 -T^2/2; 0 1 T -T; 0 0 1 0; 0 0 0 1];
 * B=[0 0;0 0; 1 0; 0 1];
 * C=[1 0 0 0; 0 0 1 0];
 * Q=[0.36 0; 0 0.05];
 * R=[100 0; 0 70]; //R=[15^2 0; 0 36];
 * Plant=ss(A,B,C,0,T,'inputname',{'v_acc_noise', 'v_acc_bias_noise'},'outputname',{'h','vAcc'},'statename',{'h','RoC','vAcc','v_acc_bias'});
 * [kalmf,L,P,M] = kalman(Plant,Q,R);
 * M
 */

/* MATLAB CODE for KF with LIDAR
 * T=1/200;
 * A=[1 T T^2/2 -T^2/2; 0 1 T -T; 0 0 1 0; 0 0 0 1];
 * B=[0 0;0 0; 1 0; 0 1];
 * C=[1 0 0 0; 0 1 0 0; 0 0 1 0];
 * Q=[0.36 0; 0 0.05];
 * R=[100 0 0; 0 80 0; 0 0 70]; //R=[15^2 0 0; 0 80 0; 0 0 36];
 * Plant=ss(A,B,C,0,T,'inputname',{'v_acc_noise', 'v_acc_bias_noise'},'outputname',{'h','RoC','vAcc'},'statename',{'h','RoC','vAcc','v_acc_bias'});
 * [kalmf,L,P,M] = kalman(Plant,Q,R);
 * M
 */

/* MATLAB CODE for KF with LIDAR and Velocity Down correction
 * T=1/200;
 * A=[1 T T^2/2 -T^2/2; 0 1 T -T; 0 0 1 0; 0 0 0 1];
 * B=[0 0;0 0; 1 0; 0 1];
 * C=[1 0 0 0; 0 1 0 0; 0 1 0 0; 0 0 1 0];
 * Q=[0.36 0; 0 0.05];
 * R=[100 0 0 0; 0 80 0 0; 0 0 80 0; 0 0 0 70];
 * Plant=ss(A,B,C,0,T,'inputname',{'v_acc_noise', 'v_acc_bias_noise'},'outputname',{'h','RoC LIDAR','RoC EKF', 'vAcc'},'statename',{'h','RoC','vAcc','v_acc_bias'});
 * [kalmf,L,P,M] = kalman(Plant,Q,R);
 * M
 */

/* Private variables ---------------------------------------------------------*/
static float _inv_pressZeroLevel; /* Ground pressure */
static float _alt_k1, _alt_kpow;  /* Altitude calculation coefficients */

#ifdef configALTITUDE_KF_ACC_HP_FILTER
IIRFilterGeneric_t HPFilt_accD;
#endif

#if (configUSE_ALT_TOF != configTOF_DISABLE)
IIRFilterDerivative_t LIDAR_diff;
#endif

/* Private functions ---------------------------------------------------------*/
static float altitudeKFAccelDownCalc(axis3f_t accel, float b_az, axis3f_t angles) {
    float accD;
#if defined(configALTITUDE_KF_USE_ACC_D)
    //transform accel_z into accel_D
    accD = (constG - accel.x * sinf(angles.y) + accel.y * cosf(angles.y) * sinf(angles.x)
            + (accel.z - b_az) * cosf(angles.y) * cosf(angles.x));
#else
    accD = (constG + accel.z - b_az);
#endif

#ifdef configALTITUDE_KF_ACC_HP_FILTER
    accD = IIRFilterProcess(&HPFilt_accD, accD);
#endif

    return DEADBAND(accD, configALTITUDE_KF_ACCEL_D_DEADBAND);
}

#ifdef configALTITUDE_KF_USE_VELD_CORRECTION
static float altitudeKFVelDownCalc(axis3f_t vel, axis3f_t angles) {
    float velD;
    //transform vel into vel_D
    velD = -vel.x * sinf(angles.y) + vel.y * cosf(angles.y) * sinf(angles.x) + vel.z * cosf(angles.y) * cosf(angles.x);

    return velD;
}
#endif

/* Functions -----------------------------------------------------------------*/

void altitudeKF_init(altitudeState_t* altState, float pressGround, float tempGround) {

    /* Initialize support variables for altitude calculation */
    _inv_pressZeroLevel = 1.f / pressGround;
    _alt_k1 = tempGround / configALTITUDE_KF_CONST_TEMP_RATE;
    _alt_kpow = (configALTITUDE_KF_CONST_R * configALTITUDE_KF_CONST_TEMP_RATE) / (configALTITUDE_KF_CONST_M * constG);

    /* Initialize state */
    altState->RoC = 0;
    altState->alt = 0;
    altState->vAcc = 0;
    altState->b_vAcc = 0;

    /* Initialize downward acceleration high-pass filter */
#ifdef configALTITUDE_KF_ACC_HP_FILTER
    IIRFilterInit(&HPFilt_accD, configALTITUDE_KF_ACCEL_HP_N0, configALTITUDE_KF_ACCEL_HP_N1,
                  configALTITUDE_KF_ACCEL_HP_N2, configALTITUDE_KF_ACCEL_HP_N3, configALTITUDE_KF_ACCEL_HP_D1,
                  configALTITUDE_KF_ACCEL_HP_D2, configALTITUDE_KF_ACCEL_HP_D3);
#endif

    /* Initialize derivative calculation for vertical speed estimation */
#if (configUSE_ALT_TOF != configTOF_DISABLE)
    IIRFilterDerivativeInit(&LIDAR_diff, configALTITUDE_KF_LIDAR_DIFF_ND, configALTITUDE_KF_LIDAR_UPDATE_TIME_S * 1e3f);
#endif
    return;
}

void altitudeKF_prediction(altitudeState_t* altState) {
    /* Predict state */
    altState->alt +=
        configALTITUDE_KF_LOOP_TIME_S * altState->RoC
        + 0.5 * configALTITUDE_KF_LOOP_TIME_S * configALTITUDE_KF_LOOP_TIME_S * (altState->vAcc - altState->b_vAcc);

    altState->RoC += configALTITUDE_KF_LOOP_TIME_S * (altState->vAcc - altState->b_vAcc);

    return;
}

void altitudeKF_updateBaroAccel(altitudeState_t* altState, float press, axis3f_t accel, float b_az, axis3f_t angles) {
    /* Calculate delta measures */
    //float delta_baroAltitude =  (1.0f - powf(press *_inv_pressZeroLevel, 0.190295f)) * 44330.0f - altState->alt;
    float delta_baroAltitude = (_alt_k1 * (1 - powf((press * _inv_pressZeroLevel), _alt_kpow))) - altState->alt;
    float delta_accelDown = altitudeKFAccelDownCalc(accel, b_az, angles) + altState->vAcc;

    /* Correct with accelerometer only if measured value is within allowed range */
    if (fabsf(delta_accelDown) > configALTITUDE_KF_MAX_ACCEL_DOWN) {
        delta_accelDown = 0;
    }

    /* Apply correction */
    altState->alt += configALTITUDE_KF_M_HH * delta_baroAltitude - configALTITUDE_KF_M_HA * delta_accelDown;
    altState->RoC += configALTITUDE_KF_M_VH * delta_baroAltitude - configALTITUDE_KF_M_VA * delta_accelDown;
    altState->vAcc += configALTITUDE_KF_M_AH * delta_baroAltitude - configALTITUDE_KF_M_AA * delta_accelDown;
    altState->b_vAcc += configALTITUDE_KF_M_BH * delta_baroAltitude - configALTITUDE_KF_M_BA * delta_accelDown;

    return;
}

#if (configUSE_ALT_TOF != configTOF_DISABLE)
void altitudeKF_updateLIDAR(altitudeState_t* altState, float ToFAlt, axis3f_t angles) {
    /* Differentiate LIDAR reading to obtain vertical speed */
    IIRFilterDerivativeProcess(&LIDAR_diff, (ToFAlt * cosf(angles.y) * cosf(angles.x)));

    /* Correct with LIDAR only if measured altitude and current attitude are within allowed range */
    if ((LIDAR_diff.output < configALTITUDE_KF_MAX_LIDAR_ROC)
        && (fabsf(angles.x) <= configALTITUDE_KF_MAX_LIDAR_ROLL_PITCH)
        && (fabsf(angles.y) <= configALTITUDE_KF_MAX_LIDAR_ROLL_PITCH)) {
        float delta_LIDARRoC =
            (LIDAR_diff.output - altState->RoC) * configALTITUDE_KF_LIDAR_UPDATE_TIME_S / configALTITUDE_KF_LOOP_TIME_S;
        altState->alt += configALTITUDE_KF_M_HL * delta_LIDARRoC;
        altState->RoC += configALTITUDE_KF_M_VL * delta_LIDARRoC;
        altState->vAcc += configALTITUDE_KF_M_AL * delta_LIDARRoC;
        altState->b_vAcc += configALTITUDE_KF_M_BL * delta_LIDARRoC;
    }
}
#endif

#ifdef configALTITUDE_KF_USE_VELD_CORRECTION
void altitudeKF_updateVelD(altitudeState_t* altState, axis3f_t velocities, axis3f_t angles) {
    float delta_velD = (altitudeKFVelDownCalc(velocities, angles) + altState->RoC);
    altState->alt -= configALTITUDE_KF_M_HV * delta_velD;
    altState->RoC -= configALTITUDE_KF_M_VV * delta_velD;
    altState->vAcc -= configALTITUDE_KF_M_AV * delta_velD;
    altState->b_vAcc -= configALTITUDE_KF_M_BV * delta_velD;
}
#endif

void altitudeKF_reset(altitudeState_t* altState, float pressGround) {
    _inv_pressZeroLevel = 1.f / pressGround;

    /* Initialize state */
    altState->RoC = 0;
    altState->alt = 0;
    altState->vAcc = 0;
    altState->b_vAcc = 0;

    /* Initialize accelerometer HP filter */
#ifdef configALTITUDE_KF_ACC_HP_FILTER
    IIRFilterReset(&HPFilt_accD);
#endif

    /* Initialize LIDAR derivative filter */
#if (configUSE_ALT_TOF != configTOF_DISABLE)
    IIRFilterDerivativeReset(&LIDAR_diff);
#endif
}
