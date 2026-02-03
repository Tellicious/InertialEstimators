/* BEGIN Header */
/**
 ******************************************************************************
 * \file            AHRS_PX4_EKF.c
 * \author          Andrea Vivani
 * \brief           PX4 attitude and heading EKF estimator
 ******************************************************************************
 * \copyright
 *
 * Copyright 2024 Andrea Vivani
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

#include "AHRS_PX4_EKF.h"
#include "basicMath.h"
#include "matrix.h"
#include "numMethods.h"

/* Private variables ---------------------------------------------------------*/

static matrix_t
    AHRS_EKF_u; //(Roll=Phi, Pitch=Theta, Yaw=Psi, Xd, Yd, Zd, c_damp, b_az, incl, bconstGx, bconstGy, bconstGz) angles in rad, angular velocities in rad/s, velocities in m/s, c_damp in N*s/m
static matrix_t _A; //state matrix
static matrix_t _C; // output matrix
static matrix_t _J, _Ji, _Q;
static matrix_t _P;      //expected errors value matrix
static matrix_t _R_acc;  //accelerometer noises covariance matrix
static matrix_t _R_gyro; //gyroscope noises covariance matrix
static matrix_t _R_mag;  //magnetometer noises covariance matrix
static matrix_t _M;      //temporary matrix
static matrix_t _K;      //gain matrix

/* Support and temporary variables */
static matrix_t TMP1, TMP2, TMP3, TMP4, TMP5;

/* Functions -----------------------------------------------------------------*/
void AHRS_PX4_EKF_init() {

    /* Initialize matrices */
    matrixInit(&AHRS_EKF_u, 12, 1);
    matrixInit(&_A, 12, 12);
    matrixInit(&_C, 3, 12);
    matrixInit(&_J, 3, 3);
    matrixInit(&_Ji, 3, 3);
    matrixInit(&_Q, 12, 12);
    matrixInit(&_P, 12, 12);
    matrixInit(&_R_acc, 3, 3);
    matrixInit(&_R_gyro, 3, 3);
    matrixInit(&_R_mag, 3, 3);
    matrixInit(&_M, 3, 3);
    matrixInit(&_K, 12, 3);
    matrixInit(&TMP1, 12, 3);
    matrixInit(&TMP2, 3, 3);
    matrixInit(&TMP3, 12, 1);
    matrixInit(&TMP4, 12, 12);
    matrixInit(&TMP5, 12, 12);

    matrixIdentity(&_P);
    matrixMultScalar(&_P, 200, &_P);
    matrixIdentity(&_J);
    matrixInversed(&_J, &_Ji);

    ELEM(AHRS_EKF_u, 8, 0) = -constG;
    ELEM(AHRS_EKF_u, 9, 0) = COS(configAHRS_PX4_EKF_INCL0);
    ELEM(AHRS_EKF_u, 11, 0) = SIN(configAHRS_PX4_EKF_INCL0);
    ELEM(_Q, 0, 0) = configAHRS_PX4_EKF_R_S_NOISE * configAHRS_PX4_EKF_PRED_LOOP_TIME_S;
    ELEM(_Q, 1, 1) = configAHRS_PX4_EKF_R_S_NOISE * configAHRS_PX4_EKF_PRED_LOOP_TIME_S;
    ELEM(_Q, 2, 2) = configAHRS_PX4_EKF_R_S_NOISE * configAHRS_PX4_EKF_PRED_LOOP_TIME_S;
    ELEM(_Q, 3, 3) = configAHRS_PX4_EKF_R_A_NOISE * configAHRS_PX4_EKF_PRED_LOOP_TIME_S;
    ELEM(_Q, 4, 4) = configAHRS_PX4_EKF_R_A_NOISE * configAHRS_PX4_EKF_PRED_LOOP_TIME_S;
    ELEM(_Q, 5, 5) = configAHRS_PX4_EKF_R_A_NOISE * configAHRS_PX4_EKF_PRED_LOOP_TIME_S;
    ELEM(_Q, 6, 6) = configAHRS_PX4_EKF_A_NOISE * configAHRS_PX4_EKF_PRED_LOOP_TIME_S;
    ELEM(_Q, 7, 7) = configAHRS_PX4_EKF_A_NOISE * configAHRS_PX4_EKF_PRED_LOOP_TIME_S;
    ELEM(_Q, 8, 8) = configAHRS_PX4_EKF_A_NOISE * configAHRS_PX4_EKF_PRED_LOOP_TIME_S;
    ELEM(_Q, 9, 9) = configAHRS_PX4_EKF_M_NOISE * configAHRS_PX4_EKF_PRED_LOOP_TIME_S;
    ELEM(_Q, 10, 10) = configAHRS_PX4_EKF_M_NOISE * configAHRS_PX4_EKF_PRED_LOOP_TIME_S;
    ELEM(_Q, 11, 11) = configAHRS_PX4_EKF_M_NOISE * configAHRS_PX4_EKF_PRED_LOOP_TIME_S;
    ELEM(_R_gyro, 0, 0) = configAHRS_PX4_EKF_GYRO_NOISE / configAHRS_PX4_EKF_GYRO_LOOP_TIME_S;
    ELEM(_R_gyro, 1, 1) = configAHRS_PX4_EKF_GYRO_NOISE / configAHRS_PX4_EKF_GYRO_LOOP_TIME_S;
    ELEM(_R_gyro, 2, 2) = configAHRS_PX4_EKF_GYRO_NOISE / configAHRS_PX4_EKF_GYRO_LOOP_TIME_S;
    ELEM(_R_acc, 0, 0) = configAHRS_PX4_EKF_ACCEL_NOISE / configAHRS_PX4_EKF_ACC_LOOP_TIME_S;
    ELEM(_R_acc, 1, 1) = configAHRS_PX4_EKF_ACCEL_NOISE / configAHRS_PX4_EKF_ACC_LOOP_TIME_S;
    ELEM(_R_acc, 2, 2) = configAHRS_PX4_EKF_ACCEL_NOISE / configAHRS_PX4_EKF_ACC_LOOP_TIME_S;
    ELEM(_R_mag, 0, 0) = configAHRS_PX4_EKF_MAG_NOISE / configAHRS_PX4_EKF_MAG_LOOP_TIME_S;
    ELEM(_R_mag, 1, 1) = configAHRS_PX4_EKF_MAG_NOISE / configAHRS_PX4_EKF_MAG_LOOP_TIME_S;
    ELEM(_R_mag, 2, 2) = configAHRS_PX4_EKF_MAG_NOISE / configAHRS_PX4_EKF_MAG_LOOP_TIME_S;

    return;
}

void AHRS_PX4_EKF_prediction() {
    /*
     wx=  u(0);   % x  body angular rate
     wy=  u(1);   % y  body angular rate
     wz=  u(2);   % z  body angular rate
     
     wax=  u(3);  % x  body angular acceleration
     way=  u(4);  % y  body angular acceleration
     waz=  u(5);  % z  body angular acceleration
     
     zex=  u(6);  % x  component gravity vector
     zey=  u(7);  % y  component gravity vector
     zez=  u(8);  % z  component gravity vector
     
     mux=  u(9); % x  component magnetic field vector
     muy=  u(10); % y  component magnetic field vector
     muz=  u(11); % z  component magnetic field vector */

    float w1 = configAHRS_PX4_EKF_PRED_LOOP_TIME_S * ELEM(AHRS_EKF_u, 0, 0);
    float w2 = configAHRS_PX4_EKF_PRED_LOOP_TIME_S * ELEM(AHRS_EKF_u, 1, 0);
    float w3 = configAHRS_PX4_EKF_PRED_LOOP_TIME_S * ELEM(AHRS_EKF_u, 2, 0);
    /* A matrix */
    ELEM(_A, 0, 0) = 1.f;
    ELEM(_A, 0, 3) = configAHRS_PX4_EKF_PRED_LOOP_TIME_S;
    ELEM(_A, 1, 1) = 1.f;
    ELEM(_A, 1, 4) = configAHRS_PX4_EKF_PRED_LOOP_TIME_S;
    ELEM(_A, 2, 2) = 1.f;
    ELEM(_A, 2, 5) = configAHRS_PX4_EKF_PRED_LOOP_TIME_S;
    ELEM(_A, 3, 3) = 1.f;
    ELEM(_A, 4, 4) = 1.f;
    ELEM(_A, 5, 5) = 1.f;
    ELEM(_A, 6, 1) = -configAHRS_PX4_EKF_PRED_LOOP_TIME_S * ELEM(AHRS_EKF_u, 8, 0);
    ELEM(_A, 6, 2) = configAHRS_PX4_EKF_PRED_LOOP_TIME_S * ELEM(AHRS_EKF_u, 7, 0);
    ELEM(_A, 6, 6) = 1.f;
    ELEM(_A, 6, 7) = w3;
    ELEM(_A, 6, 8) = -w2;
    ELEM(_A, 7, 0) = configAHRS_PX4_EKF_PRED_LOOP_TIME_S * ELEM(AHRS_EKF_u, 8, 0);
    ELEM(_A, 7, 2) = -configAHRS_PX4_EKF_PRED_LOOP_TIME_S * ELEM(AHRS_EKF_u, 6, 0);
    ELEM(_A, 7, 6) = -w3;
    ELEM(_A, 7, 7) = 1.f;
    ELEM(_A, 7, 8) = w1;
    ELEM(_A, 8, 0) = -configAHRS_PX4_EKF_PRED_LOOP_TIME_S * ELEM(AHRS_EKF_u, 7, 0);
    ELEM(_A, 8, 1) = configAHRS_PX4_EKF_PRED_LOOP_TIME_S * ELEM(AHRS_EKF_u, 6, 0);
    ELEM(_A, 8, 6) = w2;
    ELEM(_A, 8, 7) = -w1;
    ELEM(_A, 8, 8) = 1.f;
    ELEM(_A, 9, 1) = -configAHRS_PX4_EKF_PRED_LOOP_TIME_S * ELEM(AHRS_EKF_u, 11, 0);
    ELEM(_A, 9, 2) = configAHRS_PX4_EKF_PRED_LOOP_TIME_S * ELEM(AHRS_EKF_u, 10, 0);
    ELEM(_A, 9, 9) = 1.f;
    ELEM(_A, 9, 10) = w3;
    ELEM(_A, 9, 11) = -w2;
    ELEM(_A, 10, 0) = configAHRS_PX4_EKF_PRED_LOOP_TIME_S * ELEM(AHRS_EKF_u, 11, 0);
    ELEM(_A, 10, 2) = -configAHRS_PX4_EKF_PRED_LOOP_TIME_S * ELEM(AHRS_EKF_u, 9, 0);
    ELEM(_A, 10, 9) = -w3;
    ELEM(_A, 10, 10) = 1.f;
    ELEM(_A, 10, 11) = w1;
    ELEM(_A, 11, 0) = -configAHRS_PX4_EKF_PRED_LOOP_TIME_S * ELEM(AHRS_EKF_u, 10, 0);
    ELEM(_A, 11, 1) = configAHRS_PX4_EKF_PRED_LOOP_TIME_S * ELEM(AHRS_EKF_u, 9, 0);
    ELEM(_A, 11, 9) = w2;
    ELEM(_A, 11, 10) = -w1;
    ELEM(_A, 11, 11) = 1.f;

    /* Predicted P matrix */
    QuadProd(&_A, &_P, &TMP4);
    matrixCopy(&TMP4, &_P);
    matrixAdd(&_Q, &_P, &_P);

    /* Compute the apriori state estimate from the previous aposteriori estimate */

    /* Body angular rates prediction */
#if (configAHRS_PX4_EKF_INERTIA_MATRIX == configAHRS_PX4_EKF_USE_COMPLETE_INERTIA_MATRIX)

    float t1 = ELEM(_J, 0, 1) * ELEM(AHRS_EKF_u, 3, 0) + ELEM(_J, 1, 1) * ELEM(AHRS_EKF_u, 4, 0)
               + ELEM(_J, 1, 2) * ELEM(AHRS_EKF_u, 5, 0);
    float t2 = ELEM(_J, 0, 0) * ELEM(AHRS_EKF_u, 3, 0) + ELEM(_J, 0, 1) * ELEM(AHRS_EKF_u, 4, 0)
               + ELEM(_J, 0, 2) * ELEM(AHRS_EKF_u, 5, 0);
    float t3 = ELEM(_J, 0, 2) * ELEM(AHRS_EKF_u, 3, 0) + ELEM(_J, 1, 2) * ELEM(AHRS_EKF_u, 4, 0)
               + ELEM(_J, 2, 2) * ELEM(AHRS_EKF_u, 5, 0);
    float tmp1 = ELEM(AHRS_EKF_u, 3, 0) * t1 - ELEM(AHRS_EKF_u, 4, 0) * t2;
    float tmp2 = ELEM(AHRS_EKF_u, 3, 0) * t3 - ELEM(AHRS_EKF_u, 5, 0) * t2;
    float tmp3 = ELEM(AHRS_EKF_u, 4, 0) * t3 - ELEM(AHRS_EKF_u, 5, 0) * t1;

    float wax = ELEM(AHRS_EKF_u, 3, 0)
                - configAHRS_PX4_EKF_PRED_LOOP_TIME_S
                      * (ELEM(_Ji, 0, 2) * tmp1 - ELEM(_Ji, 0, 1) * tmp2 + ELEM(_Ji, 0, 0) * tmp3);
    float way = ELEM(AHRS_EKF_u, 4, 0)
                - configAHRS_PX4_EKF_PRED_LOOP_TIME_S
                      * (ELEM(_Ji, 1, 2) * tmp1 - ELEM(_Ji, 1, 1) * tmp2 + ELEM(_Ji, 0, 1) * tmp3);
    float waz = ELEM(AHRS_EKF_u, 5, 0)
                - configAHRS_PX4_EKF_PRED_LOOP_TIME_S
                      * (ELEM(_Ji, 2, 2) * tmp1 - ELEM(_Ji, 1, 2) * tmp2 + ELEM(_Ji, 0, 2) * tmp3);

#elif (configAHRS_PX4_EKF_INERTIA_MATRIX == configAHRS_PX4_EKF_USE_DIAGONAL_INERTIA_MATRIX)

    float wax = ELEM(AHRS_EKF_u, 3, 0)
                + (configAHRS_PX4_EKF_PRED_LOOP_TIME_S * ELEM(AHRS_EKF_u, 4, 0) * ELEM(AHRS_EKF_u, 5, 0)
                   * (ELEM(_J, 1, 1) - ELEM(_J, 2, 2)))
                      / ELEM(_J, 0, 0);
    float way = ELEM(AHRS_EKF_u, 4, 0)
                - (configAHRS_PX4_EKF_PRED_LOOP_TIME_S * ELEM(AHRS_EKF_u, 3, 0) * ELEM(AHRS_EKF_u, 5, 0)
                   * (ELEM(_J, 0, 0) - ELEM(_J, 2, 2)))
                      / ELEM(_J, 1, 1);
    float waz = ELEM(AHRS_EKF_u, 5, 0)
                + (configAHRS_PX4_EKF_PRED_LOOP_TIME_S * ELEM(AHRS_EKF_u, 3, 0) * ELEM(AHRS_EKF_u, 4, 0)
                   * (ELEM(_J, 0, 0) - ELEM(_J, 1, 1)))
                      / ELEM(_J, 2, 2);

#else

    float wax = ELEM(AHRS_EKF_u, 3, 0);
    float way = ELEM(AHRS_EKF_u, 4, 0);
    float waz = ELEM(AHRS_EKF_u, 5, 0);

#endif

    float delta_u6 =
        configAHRS_PX4_EKF_PRED_LOOP_TIME_S
        * (ELEM(AHRS_EKF_u, 2, 0) * ELEM(AHRS_EKF_u, 7, 0) - ELEM(AHRS_EKF_u, 1, 0) * ELEM(AHRS_EKF_u, 8, 0));

    float delta_u7 =
        configAHRS_PX4_EKF_PRED_LOOP_TIME_S
        * (ELEM(AHRS_EKF_u, 0, 0) * ELEM(AHRS_EKF_u, 8, 0) - ELEM(AHRS_EKF_u, 2, 0) * ELEM(AHRS_EKF_u, 6, 0));
    float delta_u8 =
        configAHRS_PX4_EKF_PRED_LOOP_TIME_S
        * (ELEM(AHRS_EKF_u, 1, 0) * ELEM(AHRS_EKF_u, 6, 0) - ELEM(AHRS_EKF_u, 0, 0) * ELEM(AHRS_EKF_u, 7, 0));

    float delta_u9 =
        configAHRS_PX4_EKF_PRED_LOOP_TIME_S
        * (ELEM(AHRS_EKF_u, 2, 0) * ELEM(AHRS_EKF_u, 10, 0) - ELEM(AHRS_EKF_u, 1, 0) * ELEM(AHRS_EKF_u, 11, 0));
    float delta_u10 =
        configAHRS_PX4_EKF_PRED_LOOP_TIME_S
        * (ELEM(AHRS_EKF_u, 0, 0) * ELEM(AHRS_EKF_u, 11, 0) - ELEM(AHRS_EKF_u, 2, 0) * ELEM(AHRS_EKF_u, 9, 0));
    float delta_u11 =
        configAHRS_PX4_EKF_PRED_LOOP_TIME_S
        * (ELEM(AHRS_EKF_u, 1, 0) * ELEM(AHRS_EKF_u, 9, 0) - ELEM(AHRS_EKF_u, 0, 0) * ELEM(AHRS_EKF_u, 10, 0));

    /* Predict state */
    ELEM(AHRS_EKF_u, 0, 0) += configAHRS_PX4_EKF_PRED_LOOP_TIME_S * wax;
    ELEM(AHRS_EKF_u, 1, 0) += configAHRS_PX4_EKF_PRED_LOOP_TIME_S * way;
    ELEM(AHRS_EKF_u, 2, 0) += configAHRS_PX4_EKF_PRED_LOOP_TIME_S * waz;
    ELEM(AHRS_EKF_u, 3, 0) = wax;
    ELEM(AHRS_EKF_u, 4, 0) = way;
    ELEM(AHRS_EKF_u, 5, 0) = waz;
    ELEM(AHRS_EKF_u, 6, 0) += delta_u6;
    ELEM(AHRS_EKF_u, 7, 0) += delta_u7;
    ELEM(AHRS_EKF_u, 8, 0) += delta_u8;
    ELEM(AHRS_EKF_u, 9, 0) += delta_u9;
    ELEM(AHRS_EKF_u, 10, 0) += delta_u10;
    ELEM(AHRS_EKF_u, 11, 0) += delta_u11;

    return;
}

//---------------------------------Update with gyro-----------------------------------------//
void AHRS_PX4_EKF_updateGyro(axis3f_t gyro) {
    matrix_t deltaM;
    matrixInit(&deltaM, 3, 1);

    /* C matrix */
    matrixZeros(&_C);
    ELEM(_C, 0, 0) = 1.f;
    ELEM(_C, 1, 1) = 1.f;
    ELEM(_C, 2, 2) = 1.f;

    /* Delta measures */
    ELEM(deltaM, 0, 0) = gyro.x - ELEM(AHRS_EKF_u, 0, 0);
    ELEM(deltaM, 1, 0) = gyro.y - ELEM(AHRS_EKF_u, 1, 0);
    ELEM(deltaM, 2, 0) = gyro.z - ELEM(AHRS_EKF_u, 2, 0);

    /* Gain matrix K */
    /*  _M = QuadProd(_C, _P) + _R_gyro; */
    QuadProd(&_C, &_P, &_M);
    matrixAdd(&_M, &_R_gyro, &_M);
    /* _K = _P * (~_C) * (!_M); */
    matrixMult_rhsT(&_P, &_C, &TMP1); //TMP1 contains _P * (~_C)
    matrixInversed(&_M, &TMP2);       //TMP2 contains (!_M)
    matrixMult(&TMP1, &TMP2, &_K);

    /* Correct state vector */
    //u += _K * deltaM
    matrixMult(&_K, &deltaM, &TMP3);
    matrixAdd(&AHRS_EKF_u, &TMP3, &AHRS_EKF_u);

    /* Updated P matrix */
    //_P -= _K * _C * _P;
    //_P=(_P+(~_P))*0.5; //guarantees P to be symmetric
    matrixMult(&_K, &_C, &TMP4);
    matrixMult(&TMP4, &_P, &TMP5);
    matrixSub(&_P, &TMP5, &_P);

    matrixDelete(&deltaM);

    return;
}

//---------------------------------Update with accel----------------------------------------//
void AHRS_PX4_EKF_updateAccel(axis3f_t accel) {
    matrix_t deltaM;
    matrixInit(&deltaM, 3, 1);

    /* C matrix */
    matrixZeros(&_C);
    ELEM(_C, 0, 6) = 1.f;
    ELEM(_C, 1, 7) = 1.f;
    ELEM(_C, 2, 8) = 1.f;

    /* Delta measures */
    ELEM(deltaM, 0, 0) = accel.x - ELEM(AHRS_EKF_u, 6, 0);
    ELEM(deltaM, 1, 0) = accel.y - ELEM(AHRS_EKF_u, 7, 0);
    ELEM(deltaM, 2, 0) = accel.z - ELEM(AHRS_EKF_u, 8, 0);

    /* Gain matrix K */
    /*  _M = QuadProd(_C, _P) + _R_gyro; */
    QuadProd(&_C, &_P, &_M);
    matrixAdd(&_M, &_R_acc, &_M);
    /* _K = _P * (~_C) * (!_M); */
    matrixMult_rhsT(&_P, &_C, &TMP1); //TMP1 contains _P * (~_C)
    matrixInversed(&_M, &TMP2);       //TMP2 contains (!_M)
    matrixMult(&TMP1, &TMP2, &_K);

    /* Correct state vector */
    //u += _K * deltaM
    matrixMult(&_K, &deltaM, &TMP3);
    matrixAdd(&AHRS_EKF_u, &TMP3, &AHRS_EKF_u);

    /* Updated P matrix */
    //_P -= _K * _C * _P;
    //_P=(_P+(~_P))*0.5; //guarantees P to be symmetric
    matrixMult(&_K, &_C, &TMP4);
    matrixMult(&TMP4, &_P, &TMP5);
    matrixSub(&_P, &TMP5, &_P);

    matrixDelete(&deltaM);

    return;
}

//----------------------------------Update with mag-----------------------------------------//
void AHRS_PX4_EKF_updateMag(axis3f_t mag) {

    matrix_t deltaM;
    matrixInit(&deltaM, 3, 1);

    /* C matrix */
    matrixZeros(&_C);
    ELEM(_C, 0, 9) = 1.f;
    ELEM(_C, 1, 10) = 1.f;
    ELEM(_C, 2, 11) = 1.f;

    /* Delta measures */
    ELEM(deltaM, 0, 0) = mag.x - ELEM(AHRS_EKF_u, 9, 0);
    ELEM(deltaM, 1, 0) = mag.y - ELEM(AHRS_EKF_u, 10, 0);
    ELEM(deltaM, 2, 0) = mag.z - ELEM(AHRS_EKF_u, 11, 0);

    /* Gain matrix K */
    /*  _M = QuadProd(_C, _P) + _R_gyro; */
    QuadProd(&_C, &_P, &_M);
    matrixAdd(&_M, &_R_mag, &_M);
    /* _K = _P * (~_C) * (!_M); */
    matrixMult_rhsT(&_P, &_C, &TMP1); //TMP1 contains _P * (~_C)
    matrixInversed(&_M, &TMP2);       //TMP2 contains (!_M)
    matrixMult(&TMP1, &TMP2, &_K);

    /* Correct state vector */
    //u += _K * deltaM
    matrixMult(&_K, &deltaM, &TMP3);
    matrixAdd(&AHRS_EKF_u, &TMP3, &AHRS_EKF_u);

    /* Updated P matrix */
    //_P -= _K * _C * _P;
    //_P=(_P+(~_P))*0.5; //guarantees P to be symmetric
    matrixMult(&_K, &_C, &TMP4);
    matrixMult(&TMP4, &_P, &TMP5);
    matrixSub(&_P, &TMP5, &_P);

    matrixDelete(&deltaM);

    return;
}

//--------------------------------Compute euler angles--------------------------------------//
void AHRS_PX4_EKF_calculateAngles(axis3f_t* angles) {

    /* Normalize accelerometer estimate */
    float inv_norm =
        INVSQRT(ELEM(AHRS_EKF_u, 6, 0) * ELEM(AHRS_EKF_u, 6, 0) + ELEM(AHRS_EKF_u, 7, 0) * ELEM(AHRS_EKF_u, 7, 0)
                + ELEM(AHRS_EKF_u, 8, 0) * ELEM(AHRS_EKF_u, 8, 0));
    if (isnan(inv_norm) || isinf(inv_norm)) {
        inv_norm = 1.f / constG;
    }
    float ax = ELEM(AHRS_EKF_u, 6, 0) * inv_norm;
    float ay = ELEM(AHRS_EKF_u, 7, 0) * inv_norm;
    float az = ELEM(AHRS_EKF_u, 8, 0) * inv_norm;

    /* Normalize magnetometer estimate */
    inv_norm =
        INVSQRT(ELEM(AHRS_EKF_u, 9, 0) * ELEM(AHRS_EKF_u, 9, 0) + ELEM(AHRS_EKF_u, 10, 0) * ELEM(AHRS_EKF_u, 10, 0)
                + ELEM(AHRS_EKF_u, 11, 0) * ELEM(AHRS_EKF_u, 11, 0));
    if (isnan(inv_norm) || isinf(inv_norm)) {
        inv_norm = 1.f;
    }
    float mx = ELEM(AHRS_EKF_u, 9, 0) * inv_norm;
    float my = ELEM(AHRS_EKF_u, 10, 0) * inv_norm;
    float mz = ELEM(AHRS_EKF_u, 11, 0) * inv_norm;

    /* Compute Euler angles */
    angles->x = atan2f(-ay, -az);
    angles->y = asinf(ax);
    float sPhi = SIN(angles->x);
    float cPhi = COS(angles->x);
    float sTheta = ax;
    float cTheta = COS(angles->y);
    float Yh = my * cPhi - mz * sPhi;
    float Xh = mx * cTheta + (my * sPhi + mz * cPhi) * sTheta;
    angles->z = atan2f(-Yh, Xh);

    return;
}

//--------------------------------Set inclination angle--------------------------------------//
void AHRS_Attitude_PX4_EKF_setInclination(float incl_angle) {

    ELEM(AHRS_EKF_u, 9, 0) = COS(incl_angle);
    ELEM(AHRS_EKF_u, 10, 0) = 0;
    ELEM(AHRS_EKF_u, 11, 0) = SIN(incl_angle);

    return;
}

//----------------------------------Set process noises---------------------------------------//
void AHRS_Attitude_PX4_EKF_setProcessNoise(float rot_sp, float rot_acc, float acc, float mag) {

    ELEM(_Q, 0, 0) = rot_sp * configAHRS_PX4_EKF_PRED_LOOP_TIME_S;
    ELEM(_Q, 1, 1) = rot_sp * configAHRS_PX4_EKF_PRED_LOOP_TIME_S;
    ELEM(_Q, 2, 2) = rot_sp * configAHRS_PX4_EKF_PRED_LOOP_TIME_S;
    ELEM(_Q, 3, 3) = rot_acc * configAHRS_PX4_EKF_PRED_LOOP_TIME_S;
    ELEM(_Q, 4, 4) = rot_acc * configAHRS_PX4_EKF_PRED_LOOP_TIME_S;
    ELEM(_Q, 5, 5) = rot_acc * configAHRS_PX4_EKF_PRED_LOOP_TIME_S;
    ELEM(_Q, 6, 6) = acc * configAHRS_PX4_EKF_PRED_LOOP_TIME_S;
    ELEM(_Q, 7, 7) = acc * configAHRS_PX4_EKF_PRED_LOOP_TIME_S;
    ELEM(_Q, 8, 8) = acc * configAHRS_PX4_EKF_PRED_LOOP_TIME_S;
    ELEM(_Q, 9, 9) = mag * configAHRS_PX4_EKF_PRED_LOOP_TIME_S;
    ELEM(_Q, 10, 10) = mag * configAHRS_PX4_EKF_PRED_LOOP_TIME_S;
    ELEM(_Q, 11, 11) = mag * configAHRS_PX4_EKF_PRED_LOOP_TIME_S;

    return;
}

//------------------------------------Set gyro noises----------------------------------------//
void AHRS_Attitude_PX4_EKF_setGyroNoise(float g) {

    float inv_loop_time = 1.f / configAHRS_PX4_EKF_GYRO_LOOP_TIME_S;
    ELEM(_R_gyro, 0, 0) = g * inv_loop_time;
    ELEM(_R_gyro, 1, 1) = g * inv_loop_time;
    ELEM(_R_gyro, 2, 2) = g * inv_loop_time;

    return;
}

//-------------------------------------Set acc noises----------------------------------------//
void AHRS_Attitude_PX4_EKF_setAccelNoise(float a) {

    float inv_loop_time = 1.f / configAHRS_PX4_EKF_ACC_LOOP_TIME_S;
    ELEM(_R_acc, 0, 0) = a * inv_loop_time;
    ELEM(_R_acc, 1, 1) = a * inv_loop_time;
    ELEM(_R_acc, 2, 2) = a * inv_loop_time;

    return;
}

//-------------------------------------Set mag noises----------------------------------------//
void AHRS_Attitude_PX4_EKF_setMagNoise(float m) {

    float inv_loop_time = 1.f / configAHRS_PX4_EKF_MAG_LOOP_TIME_S;
    ELEM(_R_mag, 0, 0) = m * inv_loop_time;
    ELEM(_R_mag, 1, 1) = m * inv_loop_time;
    ELEM(_R_mag, 2, 2) = m * inv_loop_time;

    return;
}
