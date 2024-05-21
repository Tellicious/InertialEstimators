/* BEGIN Header */
/**
 ******************************************************************************
 * \file            IMU_EKF.c
 * \author          Andrea Vivani
 * \brief           Implementation of AV EKF for attitude estimation (excl. Yaw)
 ******************************************************************************
 * \copyright
 *
 * Copyright 2022 Andrea Vivani
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

#include "IMU_EKF.h"
#include "basicMath.h"
#include "matrix.h"
#include "numMethods.h"

/* Private variables ---------------------------------------------------------*/

static matrix_t
    IMU_EKF_u; //state vector (Roll=Phi, Pitch=Theta, Xd, Yd, Zd, c_damp, b_az) angles in rad, angular velocities in rad/s, velocities in m/s, c_damp in N*s/m
static matrix_t _A; //state matrix
static matrix_t _B; //input matrix
static matrix_t _C; //output matrix
static matrix_t _P; //expected errors value matrix
static matrix_t _W; //gyro and acc_z noises covariance matrix (gx, gy, gz, az, c_damp, b_az)
static matrix_t _R; //acc_x and acc_y noises covariance matrix
static matrix_t _M; //temporary matrix
static matrix_t _K; //gain matrix

/* Support and temporary variables */
static float _r_vxy, _r_vz, _r_vd; //velocities noise covariances
static matrix_t TMP1, TMP2, TMP3, TMP4, TMP5;

/* Functions -----------------------------------------------------------------*/

void IMU_EKF_init(axis3f_t* angles, axis3f_t* velocities) {

    /* Initialize matrices */
    matrixInit(&IMU_EKF_u, 7, 1);
    matrixInit(&_A, 7, 7);
    matrixInit(&_B, 7, 6);
    matrixInit(&_C, 2, 7);
    matrixInit(&_P, 7, 7);
    matrixInit(&_W, 6, 6);
    matrixInit(&_R, 2, 2);
    matrixInit(&_M, 2, 2);
    matrixInit(&_K, 7, 2);
    matrixInit(&TMP1, 7, 2);
    matrixInit(&TMP2, 2, 2);
    matrixInit(&TMP3, 7, 1);
    matrixInit(&TMP4, 7, 7);
    matrixInit(&TMP5, 7, 7);

    ELEM(IMU_EKF_u, 0, 0) = configIMU_EKF_PHI0;
    ELEM(IMU_EKF_u, 1, 0) = configIMU_EKF_THETA0;
    ELEM(IMU_EKF_u, 5, 0) = configIMU_EKF_C_DAMP0;
    ELEM(_W, 0, 0) = configIMU_EKF_GXY_NOISE * configIMU_EKF_LOOP_TIME_S;
    ELEM(_W, 1, 1) = configIMU_EKF_GXY_NOISE * configIMU_EKF_LOOP_TIME_S;
    ELEM(_W, 2, 2) = configIMU_EKF_GZ_NOISE * configIMU_EKF_LOOP_TIME_S;
    ELEM(_W, 3, 3) = configIMU_EKF_AZ_NOISE * configIMU_EKF_LOOP_TIME_S;
    ELEM(_W, 4, 4) = configIMU_EKF_C_DAMP_NOISE * configIMU_EKF_LOOP_TIME_S;
    ELEM(_W, 5, 5) = configIMU_EKF_B_AZ_NOISE * configIMU_EKF_LOOP_TIME_S;
    ELEM(_R, 0, 0) = configIMU_EKF_AXY_NOISE / configIMU_EKF_LOOP_TIME_S;
    ELEM(_R, 1, 1) = configIMU_EKF_AXY_NOISE / configIMU_EKF_LOOP_TIME_S;
    _r_vxy = configIMU_EKF_VXY_NOISE;
    _r_vz = configIMU_EKF_VZ_NOISE;
    _r_vd = configIMU_EKF_VD_NOISE;

    /* Set angles */
    angles->x = ELEM(IMU_EKF_u, 0, 0); //u(0,0) is roll according to IMU ref. frame
    angles->y = ELEM(IMU_EKF_u, 1, 0); //u(1,0) is pitch according to IMU ref. frame

    /* Set velocities */
    velocities->x = ELEM(IMU_EKF_u, 2, 0); //u(2,0) is speed along local x axis according to IMU ref. frame
    velocities->y = ELEM(IMU_EKF_u, 3, 0); //u(3,0) is speed along local y axis according to IMU ref. frame
    velocities->z = ELEM(IMU_EKF_u, 4, 0); //u(4,0) is speed along local z axis according to IMU ref. frame

    return;
}

/*------------------------------------Prediction--------------------------------------------*/
void IMU_EKF_prediction(float az, axis3f_t gyro) {
    float delta_u2, delta_u3;

    /* Trig functions */
    float sPhi = SIN(ELEM(IMU_EKF_u, 0, 0));
    float cPhi = COS(ELEM(IMU_EKF_u, 0, 0));
    float sTheta = SIN(ELEM(IMU_EKF_u, 1, 0));
    float cTheta = COS(ELEM(IMU_EKF_u, 1, 0));
    float inv_cTheta = 1.0f / cTheta;
    float tTheta = sTheta * inv_cTheta;
    float tmp1 = sPhi * gyro.y + cPhi * gyro.z;
    float tmp2 = cPhi * gyro.y - sPhi * gyro.z;

    /* A matrix */
    //_A.zeros(); //zeros or not?
    ELEM(_A, 0, 0) = 1.0f + configIMU_EKF_LOOP_TIME_S * tmp2 * tTheta;
    ELEM(_A, 0, 1) = configIMU_EKF_LOOP_TIME_S * tmp1 * inv_cTheta * inv_cTheta;
    ELEM(_A, 1, 0) = -configIMU_EKF_LOOP_TIME_S * tmp1;
    ELEM(_A, 1, 1) = 1.0f;
    ELEM(_A, 2, 1) = -configIMU_EKF_LOOP_TIME_S * constG * cTheta;
    ELEM(_A, 2, 2) = 1.0f - configIMU_EKF_LOOP_TIME_S * ELEM(IMU_EKF_u, 5, 0);
    ELEM(_A, 2, 3) = configIMU_EKF_LOOP_TIME_S * gyro.z;
    ELEM(_A, 2, 4) = -configIMU_EKF_LOOP_TIME_S * gyro.y;
    ELEM(_A, 2, 5) = -configIMU_EKF_LOOP_TIME_S * ELEM(IMU_EKF_u, 2, 0);
    ELEM(_A, 3, 0) = configIMU_EKF_LOOP_TIME_S * constG * cPhi * cTheta;
    ELEM(_A, 3, 1) = -configIMU_EKF_LOOP_TIME_S * constG * sPhi * sTheta;
    ELEM(_A, 3, 2) = -configIMU_EKF_LOOP_TIME_S * gyro.z;
    ELEM(_A, 3, 3) = 1.0f - configIMU_EKF_LOOP_TIME_S * ELEM(IMU_EKF_u, 5, 0);
    ELEM(_A, 3, 4) = configIMU_EKF_LOOP_TIME_S * gyro.x;
    ELEM(_A, 3, 5) = -configIMU_EKF_LOOP_TIME_S * ELEM(IMU_EKF_u, 3, 0);
    ELEM(_A, 4, 0) = -configIMU_EKF_LOOP_TIME_S * constG * cTheta * sPhi;
    ELEM(_A, 4, 1) = -configIMU_EKF_LOOP_TIME_S * constG * cPhi * sTheta;
    ELEM(_A, 4, 4) = 1.0f;
    ELEM(_A, 4, 6) = -configIMU_EKF_LOOP_TIME_S;
    ELEM(_A, 5, 5) = 1.0f;
    ELEM(_A, 6, 6) = 1.0f;

    /* B matrix */
    //_B.zeros(); //zeros or not?
    ELEM(_B, 0, 0) = 1.0f;
    ELEM(_B, 0, 1) = sPhi * tTheta;
    ELEM(_B, 0, 2) = cPhi * tTheta;
    ELEM(_B, 1, 1) = cPhi;
    ELEM(_B, 1, 2) = -sPhi;
    ELEM(_B, 2, 1) = -ELEM(IMU_EKF_u, 4, 0);
    ELEM(_B, 2, 2) = ELEM(IMU_EKF_u, 3, 0);
    ELEM(_B, 3, 0) = ELEM(IMU_EKF_u, 4, 0);
    ELEM(_B, 3, 2) = -ELEM(IMU_EKF_u, 2, 0);
    ELEM(_B, 4, 3) = 1.0f;
    ELEM(_B, 5, 4) = 1.0f;
    ELEM(_B, 6, 5) = 1.0f;

    /* Predicted P matrix */
    //P_m = A * P_p * (~A) + B * W * (~B);
    //Q=A*B*W*(~B)*(~A); //with continuous-time A (Ad=I+A*dt), it should be Q=A*B*W*(~B)*(~A)*T_samp but T_samp is already included in _W
    //Q=B*W*(~B); //it should be Q=B*W*(~B)*T_samp but T_samp is already included in _W
    //_P = QuadProd(_A, _P) + QuadProd(_B, _W);
    QuadProd(&_A, &_P, &TMP4);
    matrixCopy(&TMP4, &_P);
    QuadProd(&_B, &_W, &TMP4);
    matrixAdd(&TMP4, &_P, &_P);

#ifdef configIMU_EKF_CORRECT_ACCEL_OFFSET
    az += ((gyro.x * gyro.x + gyro.y * gyro.y) * configIMU_EKF_ACCEL_OFFSET_Z)
          - (gyro.x * gyro.y * configIMU_EKF_ACCEL_OFFSET_X) - (gyro.y * gyro.z * configIMU_EKF_ACCEL_OFFSET_Y);
#endif

    /* Predict state */
    ELEM(IMU_EKF_u, 0, 0) += configIMU_EKF_LOOP_TIME_S * (gyro.x + tmp1 * tTheta);
    ELEM(IMU_EKF_u, 1, 0) += configIMU_EKF_LOOP_TIME_S * tmp2;
    delta_u2 = configIMU_EKF_LOOP_TIME_S
               * (ELEM(IMU_EKF_u, 3, 0) * gyro.z - ELEM(IMU_EKF_u, 4, 0) * gyro.y
                  - ELEM(IMU_EKF_u, 5, 0) * ELEM(IMU_EKF_u, 2, 0) - constG * sTheta);
    delta_u3 = configIMU_EKF_LOOP_TIME_S
               * (ELEM(IMU_EKF_u, 4, 0) * gyro.x - ELEM(IMU_EKF_u, 2, 0) * gyro.z
                  - ELEM(IMU_EKF_u, 5, 0) * ELEM(IMU_EKF_u, 3, 0) + constG * sPhi * cTheta);
    ELEM(IMU_EKF_u, 2, 0) += delta_u2;
    ELEM(IMU_EKF_u, 3, 0) += delta_u3;
    ELEM(IMU_EKF_u, 4, 0) += configIMU_EKF_LOOP_TIME_S * (az - ELEM(IMU_EKF_u, 6, 0) + constG * cPhi * cTheta);
    // ELEM(IMU_EKF_u, 5,0) += 0;
    // ELEM(IMU_EKF_u, 6,0) += 0;

    return;
}

/*-----------------------------Update with accel & gyro------------------------------------*/
void IMU_EKF_updateAccelGyro(axis3f_t* angles, axis3f_t* velocities, axis3f_t accel, axis3f_t gyro) {
    matrix_t deltaM;
    matrixInit(&deltaM, 2, 1);

#ifdef configIMU_EKF_CORRECT_ACCEL_OFFSET
    accel.x += ((gyro.y * gyro.y + gyro.z * gyro.z) * configIMU_EKF_ACCEL_OFFSET_X)
               - (gyro.x * gyro.y * configIMU_EKF_ACCEL_OFFSET_Y) - (gyro.x * gyro.z * configIMU_EKF_ACCEL_OFFSET_Z);
    accel.y += ((gyro.x * gyro.x + gyro.z * gyro.z) * configIMU_EKF_ACCEL_OFFSET_Y)
               - (gyro.x * gyro.y * configIMU_EKF_ACCEL_OFFSET_X) - (gyro.y * gyro.z * configIMU_EKF_ACCEL_OFFSET_Z);
#endif

    /* C matrix */
    //_C.zeros(); //zeros or not?
    ELEM(_C, 0, 2) = -ELEM(IMU_EKF_u, 5, 0);
    ELEM(_C, 0, 3) = gyro.z;
    ELEM(_C, 0, 4) = -gyro.y;
    ELEM(_C, 0, 5) = -ELEM(IMU_EKF_u, 2, 0);
    ELEM(_C, 1, 2) = -gyro.z;
    ELEM(_C, 1, 3) = -ELEM(IMU_EKF_u, 5, 0);
    ELEM(_C, 1, 4) = gyro.x;
    ELEM(_C, 1, 5) = -ELEM(IMU_EKF_u, 3, 0);

    /* Delta measures */
    ELEM(deltaM, 0, 0) = accel.x + ELEM(IMU_EKF_u, 4, 0) * gyro.y - ELEM(IMU_EKF_u, 3, 0) * gyro.z
                         + ELEM(IMU_EKF_u, 5, 0) * ELEM(IMU_EKF_u, 2, 0);
    ELEM(deltaM, 1, 0) = accel.y + ELEM(IMU_EKF_u, 2, 0) * gyro.z - ELEM(IMU_EKF_u, 4, 0) * gyro.x
                         + ELEM(IMU_EKF_u, 5, 0) * ELEM(IMU_EKF_u, 3, 0);

    /* Gain matrix K */
    //_M = QuadProd(_C, _P) + _R;
    QuadProd(&_C, &_P, &_M);
    matrixAdd(&_M, &_R, &_M);
    //_K = _P * (~_C) * (!_M);
    matrixMult_rhsT(&_P, &_C, &TMP1); //TMP1 contains _P * (~_C)
    matrixInversed(&_M, &TMP2);       //TMP2 contains (!_M)
    matrixMult(&TMP1, &TMP2, &_K);

    /* Correct state vector */
    //u += _K * deltaM
    matrixMult(&_K, &deltaM, &TMP3);
    matrixAdd(&IMU_EKF_u, &TMP3, &IMU_EKF_u);

    /* Updated P matrix */
    //_P -= _K * _C * _P;
    //_P=(_P+(~_P))*0.5; //guarantees P to be symmetric
    matrixMult(&_K, &_C, &TMP4);
    matrixMult(&TMP4, &_P, &TMP5);
    matrixSub(&_P, &TMP5, &_P);

    /* Set angles */
    angles->x = matrixGet(&IMU_EKF_u, 0, 0); //u(0,0) is roll according to IMU ref. frame
    angles->y = matrixGet(&IMU_EKF_u, 1, 0); //u(1,0) is pitch according to IMU ref. frame

    /* Set velocities */
    velocities->x = ELEM(IMU_EKF_u, 2, 0); //u(2,0) is speed along local x axis according to IMU ref. frame
    velocities->y = ELEM(IMU_EKF_u, 3, 0); //u(3,0) is speed along local y axis according to IMU ref. frame
    velocities->z = ELEM(IMU_EKF_u, 4, 0); //u(4,0) is speed along local z axis according to IMU ref. frame

    matrixDelete(&deltaM);
    return;
}

/*------------------------Update with velocity along x,y local-----------------------------*/
void IMU_EKF_updateVelXY(axis3f_t* angles, axis3f_t* velocities, float vx, float vy, float dt_s) {
    matrix_t deltaM, R_tmp, C_tmp;

    matrixInit(&deltaM, 2, 1);
    matrixInit(&R_tmp, 2, 2);
    matrixInit(&C_tmp, 2, 7);

    /* R matrix */
    ELEM(R_tmp, 0, 0) = _r_vxy / dt_s;
    ELEM(R_tmp, 1, 1) = ELEM(R_tmp, 0, 0);

    /* C matrix */
    ELEM(C_tmp, 0, 2) = 1.f;
    ELEM(C_tmp, 1, 3) = 1.f;

    /* Delta measures */
    ELEM(deltaM, 0, 0) = vx - ELEM(IMU_EKF_u, 2, 0);
    ELEM(deltaM, 1, 0) = vy - ELEM(IMU_EKF_u, 3, 0);

    /* Gain matrix K */
    //_M = QuadProd(C_tmp,_P) + R_tmp;
    QuadProd(&C_tmp, &_P, &_M);
    matrixAdd(&_M, &R_tmp, &_M);
    //_K = _P * (~C_tmp) * (!_M);
    matrixMult_rhsT(&_P, &C_tmp, &TMP1); //TMP1 contains _P * (~_C)
    matrixInversed(&_M, &TMP2);          //TMP2 contains (!_M)
    matrixMult(&TMP1, &TMP2, &_K);

    /* Correct state vector */
    //u += _K * deltaM;
    matrixMult(&_K, &deltaM, &TMP3);
    matrixAdd(&IMU_EKF_u, &TMP3, &IMU_EKF_u);

    /* Updated P matrix */
    //_P -= _K * C_tmp * _P;
    //_P=(_P+(~_P))*0.5; //guarantees P to be symmetric
    matrixMult(&_K, &C_tmp, &TMP4);
    matrixMult(&TMP4, &_P, &TMP5);
    matrixSub(&_P, &TMP5, &_P);

    /* Set angles */
    angles->x = matrixGet(&IMU_EKF_u, 0, 0); //u(0,0) is roll according to IMU ref. frame
    angles->y = matrixGet(&IMU_EKF_u, 1, 0); //u(1,0) is pitch according to IMU ref. frame

    /* Set velocities */
    velocities->x = matrixGet(&IMU_EKF_u, 2, 0); //u(2,0) is speed along local x axis according to IMU ref. frame
    velocities->y = matrixGet(&IMU_EKF_u, 3, 0); //u(3,0) is speed along local x axis according to IMU ref. frame
    velocities->z = matrixGet(&IMU_EKF_u, 4, 0); //u(4,0) is speed along local x axis according to IMU ref. frame

    matrixDelete(&deltaM);
    matrixDelete(&R_tmp);
    matrixDelete(&C_tmp);

    return;
}

/*-------------------------Update with velocity along z local-------------------------------*/
void IMU_EKF_updateVelZ(axis3f_t* angles, axis3f_t* velocities, float vz, float dt_s) {
    matrix_t C_tmp, K;

    matrixInit(&C_tmp, 1, 7);
    matrixInit(&K, 7, 1);

    /* C matrix */
    matrixSet(&C_tmp, 0, 4, 1.f);

    /* Delta measures */
    float deltaM = vz - matrixGet(&IMU_EKF_u, 4, 0);

    /* Gain matrix K (change K with _K in correction and update of P) */
    /*_M = QuadProd(C_tmp,_P) + (_r_vz / dt_s);
   _K = _P * (~C_tmp) * (!_M);*/

    /* Faster Gain matrix K */
    float inv_m = 1.f / (ELEM(_P, 4, 4) + (_r_vz / dt_s));
    ELEM(K, 0, 0) = ELEM(_P, 0, 4) * inv_m;
    ELEM(K, 1, 0) = ELEM(_P, 1, 4) * inv_m;
    ELEM(K, 2, 0) = ELEM(_P, 2, 4) * inv_m;
    ELEM(K, 3, 0) = ELEM(_P, 3, 4) * inv_m;
    ELEM(K, 4, 0) = ELEM(_P, 4, 4) * inv_m;
    ELEM(K, 5, 0) = ELEM(_P, 5, 4) * inv_m;
    ELEM(K, 6, 0) = ELEM(_P, 6, 4) * inv_m;

    /* Correct state vector */
    //u += K * deltaM;
    matrixMultScalar(&K, deltaM, &TMP3);
    matrixAdd(&IMU_EKF_u, &TMP3, &IMU_EKF_u);

    /* Updated P matrix */
    //_P -= K * C_tmp * _P;
    //_P=(_P+(~_P))*0.5; //guarantees P to be symmetric
    matrixMult(&K, &C_tmp, &TMP4);
    matrixMult(&TMP4, &_P, &TMP5);
    matrixSub(&_P, &TMP5, &_P);

    /* Set angles */
    angles->x = matrixGet(&IMU_EKF_u, 0, 0); //u(0,0) is roll according to IMU ref. frame
    angles->y = matrixGet(&IMU_EKF_u, 1, 0); //u(1,0) is pitch according to IMU ref. frame

    /* Set velocities */
    velocities->x = matrixGet(&IMU_EKF_u, 2, 0); //u(2,0) is speed along local x axis according to IMU ref. frame
    velocities->y = matrixGet(&IMU_EKF_u, 3, 0); //u(3,0) is speed along local x axis according to IMU ref. frame
    velocities->z = matrixGet(&IMU_EKF_u, 4, 0); //u(4,0) is speed along local x axis according to IMU ref. frame

    matrixDelete(&C_tmp);
    matrixDelete(&K);

    return;
}

/*-------------------------Update with velocity along d global-------------------------------*/
void IMU_EKF_updateVelD(axis3f_t* angles, axis3f_t* velocities, float vD, float dt_s) {
    matrix_t C_tmp, K, M;

    matrixInit(&C_tmp, 1, 7);
    matrixInit(&K, 7, 1);
    matrixInit(&M, 1, 1);

    /* Trig functions */
    float sPhi = SIN(ELEM(IMU_EKF_u, 0, 0));
    float cPhi = COS(ELEM(IMU_EKF_u, 0, 0));
    float sTheta = SIN(ELEM(IMU_EKF_u, 1, 0));
    float cTheta = COS(ELEM(IMU_EKF_u, 1, 0));

    /* C matrix */
    ELEM(C_tmp, 0, 0) = ELEM(IMU_EKF_u, 3, 0) * cPhi * cTheta - ELEM(IMU_EKF_u, 4, 0) * cTheta * sPhi;
    ELEM(C_tmp, 0, 1) =
        -ELEM(IMU_EKF_u, 2, 0) * cTheta - ELEM(IMU_EKF_u, 4, 0) * cPhi * sTheta - ELEM(IMU_EKF_u, 3, 0) * sPhi * sTheta;
    ELEM(C_tmp, 0, 2) = -sTheta;
    ELEM(C_tmp, 0, 3) = cTheta * sPhi;
    ELEM(C_tmp, 0, 4) = cPhi * cTheta;
    ELEM(C_tmp, 0, 5) = 0;
    ELEM(C_tmp, 0, 6) = 0;

    /* Delta measures */
    float deltaM = -vD
                   - (ELEM(IMU_EKF_u, 4, 0) * cPhi * cTheta - ELEM(IMU_EKF_u, 2, 0) * sTheta
                      + ELEM(IMU_EKF_u, 3, 0) * cTheta * sPhi);

    /* Gain matrix K */
    //_M = QuadProd(C_tmp,_P) + (_r_vd / dt_s);
    QuadProd(&C_tmp, &_P, &M);
    ELEM(M, 0, 0) += (_r_vd / dt_s);
    //K = _P * (~C_tmp) * (!_M);
    matrixMult_rhsT(&_P, &C_tmp, &K);
    matrixMultScalar(&K, 1.0f / ELEM(M, 0, 0), &K);

    /* Correct state vector */
    //u += K * deltaM;
    matrixMultScalar(&K, deltaM, &TMP3);
    matrixAdd(&IMU_EKF_u, &TMP3, &IMU_EKF_u);

    /* Updated P matrix */
    //_P -= _K * C_tmp * _P;
    //_P=(_P+(~_P))*0.5; //guarantees P to be symmetric
    matrixMult(&K, &C_tmp, &TMP4);
    matrixMult(&TMP4, &_P, &TMP5);
    matrixSub(&_P, &TMP5, &_P);

    /* Set angles */
    angles->x = matrixGet(&IMU_EKF_u, 0, 0); //u(0,0) is roll according to IMU ref. frame
    angles->y = matrixGet(&IMU_EKF_u, 1, 0); //u(1,0) is pitch according to IMU ref. frame

    /* Set velocities */
    velocities->x = matrixGet(&IMU_EKF_u, 2, 0); //u(2,0) is speed along local x axis according to IMU ref. frame
    velocities->y = matrixGet(&IMU_EKF_u, 3, 0); //u(3,0) is speed along local x axis according to IMU ref. frame
    velocities->z = matrixGet(&IMU_EKF_u, 4, 0); //u(4,0) is speed along local x axis according to IMU ref. frame

    matrixDelete(&C_tmp);
    matrixDelete(&K);
    matrixDelete(&M);

    return;
}

/*----------------------------------Starting values----------------------------------------*/
void IMU_EKF_reset(axis3f_t* angles, axis3f_t* velocities, float phi0, float theta0) {
    matrixSet(&IMU_EKF_u, 0, 0, phi0);
    matrixSet(&IMU_EKF_u, 1, 0, theta0);
    matrixSet(&IMU_EKF_u, 2, 0, 0);
    matrixSet(&IMU_EKF_u, 3, 0, 0);
    matrixSet(&IMU_EKF_u, 4, 0, 0);
    //matrixSet(&IMU_EKF_u, 5,0, c_damp_0);

    /* Set angles */
    angles->x = matrixGet(&IMU_EKF_u, 0, 0); //u(0,0) is roll according to IMU ref. frame
    angles->y = matrixGet(&IMU_EKF_u, 1, 0); //u(1,0) is pitch according to IMU ref. frame

    /* Set velocities */
    velocities->x = matrixGet(&IMU_EKF_u, 2, 0); //u(2,0) is speed along local x axis according to IMU ref. frame
    velocities->y = matrixGet(&IMU_EKF_u, 3, 0); //u(3,0) is speed along local x axis according to IMU ref. frame
    velocities->z = matrixGet(&IMU_EKF_u, 4, 0); //u(4,0) is speed along local x axis according to IMU ref. frame

    return;
}

/*------------------------------------Input noises-----------------------------------------*/
void IMU_EKF_setInputNoises(float az, float gxy, float gz, float c_damp, float b_az) {
    matrixSet(&_W, 0, 0, gxy * configIMU_EKF_LOOP_TIME_S);
    matrixSet(&_W, 1, 1, gxy * configIMU_EKF_LOOP_TIME_S);
    matrixSet(&_W, 2, 2, gz * configIMU_EKF_LOOP_TIME_S);
    matrixSet(&_W, 3, 3, az * configIMU_EKF_LOOP_TIME_S);
    matrixSet(&_W, 4, 4, c_damp * configIMU_EKF_LOOP_TIME_S);
    matrixSet(&_W, 5, 5, b_az * configIMU_EKF_LOOP_TIME_S);

    return;
}

/*---------------------------------Output accel noises-------------------------------------*/
void IMU_EKF_setAccelNoise(float axy) {
    float inv_loop_time_s = 1.0f / configIMU_EKF_LOOP_TIME_S;
    matrixSet(&_R, 0, 0, axy * inv_loop_time_s);
    matrixSet(&_R, 1, 1, axy * inv_loop_time_s);

    return;
}

/*------------------------Output velocity noises along x,y local----------------------------*/
void IMU_EKF_setVelXYNoise(float vxy) {
    _r_vxy = vxy;

    return;
}

/*-------------------------Output velocity noises along z local-----------------------------*/
void IMU_EKF_setVelZNoise(float vz) {
    _r_vz = vz;

    return;
}

/*------------------------Output velocity noises along d global-----------------------------*/
void IMU_EKF_setVelDNoise(float vD) {
    _r_vd = vD;

    return;
}

/*----------------------------Get values from state vector----------------------------------*/
float IMU_EKF_getStateValue(uint8_t idx) { return matrixGet(&IMU_EKF_u, idx, 0); }
