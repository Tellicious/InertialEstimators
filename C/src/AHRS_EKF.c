/* BEGIN Header */
/**
 ******************************************************************************
 * \file            AHRS_EKF.c
 * \author          Andrea Vivani
 * \brief           Implementation of AV EKF for attitude and heading estimation
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

#include "AHRS_EKF.h"
#include "basicMath.h"
#include "math.h"
#include "matrix.h"
#include "numMethods.h"

/* Private variables ---------------------------------------------------------*/

static matrix_t
    AHRS_EKF_u; //(Roll=Phi, Pitch=Theta, Yaw=Psi, Xd, Yd, Zd, c_damp, b_az, incl, bconstGx, bconstGy, bconstGz) angles in rad, angular velocities in rad/s, velocities in m/s, c_damp in N*s/m
static matrix_t _A;     //state matrix
static matrix_t _B;     //input matrix
static matrix_t _C_acc; //acceleration output matrix
static matrix_t _C_mag; //magnetic field output matrix
static matrix_t _P;     //expected errors value matrix
static matrix_t
    _W; //gyro and acc_z noises covariance matrix (g_x, g_y, g_z, a_z, c_damp, b_az, incl, bconstGx, bconstGy, bconstGz)
static matrix_t _R_acc; //acc_x and acc_y noises covariance matrix
static matrix_t _R_mag; //magnetometer noises covariance matrix
static matrix_t _M;     //temporary matrix
static matrix_t _K;     //gain matrix

/* Support and temporary variables */
static float _r_vxy, _r_vz, _r_vne, _r_vd; //velocities noise covariances
static matrix_t TMP1, TMP2, TMP3, TMP4, TMP5;

/* Functions -----------------------------------------------------------------*/
void AHRS_EKF_init(axis3f_t* angles, axis3f_t* velocities) {

    /* Initialize matrices */
    matrixInit(&AHRS_EKF_u, 12, 1);
    matrixInit(&_A, 12, 12);
    matrixInit(&_B, 12, 10);
    matrixInit(&_C_acc, 2, 12);
    matrixInit(&_C_mag, 3, 12);
    matrixInit(&_P, 12, 12);
    matrixInit(&_W, 10, 10);
    matrixInit(&_R_acc, 2, 2);
    matrixInit(&_R_mag, 3, 3);
    matrixInit(&_M, 2, 2);
    matrixInit(&_K, 12, 2);
    matrixInit(&TMP1, 12, 2);
    matrixInit(&TMP2, 2, 2);
    matrixInit(&TMP3, 12, 1);
    matrixInit(&TMP4, 12, 12);
    matrixInit(&TMP5, 12, 12);

    ELEM(AHRS_EKF_u, 0, 0) = configAHRS_EKF_PHI0;
    ELEM(AHRS_EKF_u, 1, 0) = configAHRS_EKF_THETA0;
    ELEM(AHRS_EKF_u, 2, 0) = configAHRS_EKF_PSI0;
    ELEM(AHRS_EKF_u, 6, 0) = configAHRS_EKF_C_DAMP0;
    ELEM(AHRS_EKF_u, 8, 0) = configAHRS_EKF_INCL0;

    ELEM(_W, 0, 0) = configAHRS_EKF_GXY_NOISE * configAHRS_EKF_LOOP_TIME_S;
    ELEM(_W, 1, 1) = configAHRS_EKF_GXY_NOISE * configAHRS_EKF_LOOP_TIME_S;
    ELEM(_W, 2, 2) = configAHRS_EKF_GZ_NOISE * configAHRS_EKF_LOOP_TIME_S;
    ELEM(_W, 3, 3) = configAHRS_EKF_AZ_NOISE * configAHRS_EKF_LOOP_TIME_S;
    ELEM(_W, 4, 4) = configAHRS_EKF_C_DAMP_NOISE * configAHRS_EKF_LOOP_TIME_S;
    ELEM(_W, 5, 5) = configAHRS_EKF_B_AZ_NOISE * configAHRS_EKF_LOOP_TIME_S;
    ELEM(_W, 6, 6) = configAHRS_EKF_INCL_NOISE * configAHRS_EKF_LOOP_TIME_S;
    ELEM(_W, 7, 7) = configAHRS_EKF_B_G_NOISE * configAHRS_EKF_LOOP_TIME_S;
    ELEM(_W, 8, 8) = configAHRS_EKF_B_G_NOISE * configAHRS_EKF_LOOP_TIME_S;
    ELEM(_W, 9, 9) = configAHRS_EKF_B_G_NOISE * configAHRS_EKF_LOOP_TIME_S;
    ELEM(_R_acc, 0, 0) = configAHRS_EKF_AXY_NOISE / configAHRS_EKF_LOOP_TIME_S;
    ELEM(_R_acc, 1, 1) = configAHRS_EKF_AXY_NOISE / configAHRS_EKF_LOOP_TIME_S;
    ELEM(_R_mag, 0, 0) = configAHRS_EKF_M_NOISE / configAHRS_EKF_MAG_LOOP_TIME_S;
    ELEM(_R_mag, 1, 1) = configAHRS_EKF_M_NOISE / configAHRS_EKF_MAG_LOOP_TIME_S;
    ELEM(_R_mag, 2, 2) = configAHRS_EKF_M_NOISE / configAHRS_EKF_MAG_LOOP_TIME_S;

    _r_vxy = configAHRS_EKF_VXY_NOISE;
    _r_vz = configAHRS_EKF_VZ_NOISE;
    _r_vne = configAHRS_EKF_VNE_NOISE;
    _r_vd = configAHRS_EKF_VD_NOISE;

    /* Set angles */
    angles->x = ELEM(AHRS_EKF_u, 0, 0); //u(0,0) is roll according to IMU ref. frame
    angles->y = ELEM(AHRS_EKF_u, 1, 0); //u(1,0) is pitch according to IMU ref. frame
    angles->z = ELEM(AHRS_EKF_u, 2, 0); //u(2,0) is yaw according to IMU ref. frame

    /* Set velocities */
    velocities->x = ELEM(AHRS_EKF_u, 3, 0); //u(3,0) is speed along local x axis according to IMU ref. frame
    velocities->y = ELEM(AHRS_EKF_u, 4, 0); //u(4,0) is speed along local y axis according to IMU ref. frame
    velocities->z = ELEM(AHRS_EKF_u, 5, 0); //u(5,0) is speed along local z axis according to IMU ref. frame

    return;
}

/*------------------------------------Prediction--------------------------------------------*/
void AHRS_EKF_prediction(float az, axis3f_t gyro) {

    /* Trig functions */
    float sPhi = sinf(ELEM(AHRS_EKF_u, 0, 0));
    float cPhi = cosf(ELEM(AHRS_EKF_u, 0, 0));
    float sTheta = sinf(ELEM(AHRS_EKF_u, 1, 0));
    float cTheta = cosf(ELEM(AHRS_EKF_u, 1, 0));
    float inv_cTheta = 1.0f / cTheta;
    float tTheta = sTheta * inv_cTheta;
    float pr = gyro.x - ELEM(AHRS_EKF_u, 9, 0);
    float qr = gyro.y - ELEM(AHRS_EKF_u, 10, 0);
    float rr = gyro.z - ELEM(AHRS_EKF_u, 11, 0);
    float tmp1 = sPhi * qr + cPhi * rr;
    float tmp2 = cPhi * qr - sPhi * rr;

    /* A matrix */
    //_A.zeros(); //zeros or not?
    ELEM(_A, 0, 0) = 1.0f + configAHRS_EKF_LOOP_TIME_S * tmp2 * tTheta;
    ELEM(_A, 0, 1) = configAHRS_EKF_LOOP_TIME_S * tmp1 * inv_cTheta * inv_cTheta;
    ELEM(_A, 0, 9) = -configAHRS_EKF_LOOP_TIME_S;
    ELEM(_A, 0, 10) = -configAHRS_EKF_LOOP_TIME_S * sPhi * tTheta;
    ELEM(_A, 0, 11) = -configAHRS_EKF_LOOP_TIME_S * cPhi * tTheta;
    ELEM(_A, 1, 0) = -configAHRS_EKF_LOOP_TIME_S * tmp1;
    ELEM(_A, 1, 1) = 1.0f;
    ELEM(_A, 1, 10) = -configAHRS_EKF_LOOP_TIME_S * cPhi;
    ELEM(_A, 1, 11) = configAHRS_EKF_LOOP_TIME_S * sPhi;
    ELEM(_A, 2, 0) = configAHRS_EKF_LOOP_TIME_S * tmp2 * inv_cTheta;
    ELEM(_A, 2, 1) = configAHRS_EKF_LOOP_TIME_S * tmp1 * tTheta * inv_cTheta;
    ELEM(_A, 2, 2) = 1.0f;
    ELEM(_A, 2, 10) = -configAHRS_EKF_LOOP_TIME_S * sPhi * inv_cTheta;
    ELEM(_A, 2, 11) = -configAHRS_EKF_LOOP_TIME_S * cPhi * inv_cTheta;
    ELEM(_A, 3, 1) = -configAHRS_EKF_LOOP_TIME_S * constG * cTheta;
    ELEM(_A, 3, 3) = 1.0f - configAHRS_EKF_LOOP_TIME_S * ELEM(AHRS_EKF_u, 6, 0);
    ELEM(_A, 3, 4) = configAHRS_EKF_LOOP_TIME_S * rr;
    ELEM(_A, 3, 5) = -configAHRS_EKF_LOOP_TIME_S * qr;
    ELEM(_A, 3, 6) = -configAHRS_EKF_LOOP_TIME_S * ELEM(AHRS_EKF_u, 3, 0);
    ELEM(_A, 4, 0) = configAHRS_EKF_LOOP_TIME_S * constG * cPhi * cTheta;
    ELEM(_A, 4, 1) = -configAHRS_EKF_LOOP_TIME_S * constG * sPhi * sTheta;
    ELEM(_A, 4, 3) = -configAHRS_EKF_LOOP_TIME_S * rr;
    ELEM(_A, 4, 4) = 1.0f - configAHRS_EKF_LOOP_TIME_S * ELEM(AHRS_EKF_u, 6, 0);
    ELEM(_A, 4, 5) = configAHRS_EKF_LOOP_TIME_S * pr;
    ELEM(_A, 4, 6) = -configAHRS_EKF_LOOP_TIME_S * ELEM(AHRS_EKF_u, 4, 0);
    ELEM(_A, 5, 0) = configAHRS_EKF_LOOP_TIME_S * (ELEM(AHRS_EKF_u, 7, 0) - constG) * cTheta * sPhi;
    ELEM(_A, 5, 1) = configAHRS_EKF_LOOP_TIME_S * (ELEM(AHRS_EKF_u, 7, 0) - constG) * cPhi * sTheta;
    ELEM(_A, 5, 5) = 1.0f;
    ELEM(_A, 5, 7) = -configAHRS_EKF_LOOP_TIME_S * cPhi * cTheta;
    ELEM(_A, 6, 6) = 1.0f;
    ELEM(_A, 7, 7) = 1.0f;
    ELEM(_A, 8, 8) = 1.0f;
    ELEM(_A, 9, 9) = 1.0f;
    ELEM(_A, 10, 10) = 1.0f;
    ELEM(_A, 11, 11) = 1.0f;

    /* B matrix */
    //_B.zeros(); //zeros or not?
    ELEM(_B, 0, 0) = 1.0f;
    ELEM(_B, 0, 1) = sPhi * tTheta;
    ELEM(_B, 0, 2) = cPhi * tTheta;
    ELEM(_B, 1, 1) = cPhi;
    ELEM(_B, 1, 2) = -sPhi;
    ELEM(_B, 2, 1) = sPhi * inv_cTheta;
    ELEM(_B, 2, 2) = cPhi * inv_cTheta;
    ELEM(_B, 3, 1) = -ELEM(AHRS_EKF_u, 5, 0);
    ELEM(_B, 3, 2) = ELEM(AHRS_EKF_u, 4, 0);
    ELEM(_B, 4, 0) = ELEM(AHRS_EKF_u, 5, 0);
    ELEM(_B, 4, 2) = -ELEM(AHRS_EKF_u, 3, 0);
    ELEM(_B, 5, 3) = 1.0f;
    ELEM(_B, 6, 4) = 1.0f;
    ELEM(_B, 7, 5) = 1.0f;
    ELEM(_B, 8, 6) = 1.0f;
    ELEM(_B, 9, 7) = 1.0f;
    ELEM(_B, 10, 8) = 1.0f;
    ELEM(_B, 11, 9) = 1.0f;

    /* Predicted P matrix */
    //P_m = A * P_p * (~A) + B * W * (~B);
    //Q=A*B*W*(~B)*(~A); //with continuous-time A (Ad=I+A*dt), it should be Q=A*B*W*(~B)*(~A)*T_samp but T_samp is already included in _W
    //Q=B*W*(~B); //it should be Q=B*W*(~B)*T_samp but T_samp is already included in _W
    //_P = QuadProd(_A, _P) + QuadProd(_B, _W);
    QuadProd(&_A, &_P, &TMP4);
    matrixCopy(&TMP4, &_P);
    QuadProd(&_B, &_W, &TMP4);
    matrixAdd(&TMP4, &_P, &_P);

#ifdef configAHRS_EKF_CORRECT_ACCEL_OFFSET
    az += ((pr * pr + qr * qr) * configAHRS_EKF_ACCEL_OFFSET_Z) - (pr * qr * configAHRS_EKF_ACCEL_OFFSET_X)
          - (qr * rr * configAHRS_EKF_ACCEL_OFFSET_Y);
#endif

    /* Predict state */
    ELEM(AHRS_EKF_u, 0, 0) += configAHRS_EKF_LOOP_TIME_S * (pr + tmp1 * tTheta);
    ELEM(AHRS_EKF_u, 1, 0) += configAHRS_EKF_LOOP_TIME_S * tmp2;
    ELEM(AHRS_EKF_u, 2, 0) += configAHRS_EKF_LOOP_TIME_S * tmp1 * inv_cTheta;
    float delta_u3 = configAHRS_EKF_LOOP_TIME_S
                     * (ELEM(AHRS_EKF_u, 4, 0) * rr - ELEM(AHRS_EKF_u, 5, 0) * qr
                        - ELEM(AHRS_EKF_u, 6, 0) * ELEM(AHRS_EKF_u, 3, 0) - constG * sTheta);
    float delta_u4 = configAHRS_EKF_LOOP_TIME_S
                     * (ELEM(AHRS_EKF_u, 5, 0) * pr - ELEM(AHRS_EKF_u, 3, 0) * rr
                        - ELEM(AHRS_EKF_u, 6, 0) * ELEM(AHRS_EKF_u, 4, 0) + constG * sPhi * cTheta);
    ELEM(AHRS_EKF_u, 3, 0) += delta_u3;
    ELEM(AHRS_EKF_u, 4, 0) += delta_u4;
    ELEM(AHRS_EKF_u, 5, 0) += configAHRS_EKF_LOOP_TIME_S * (az + (constG - ELEM(AHRS_EKF_u, 7, 0)) * cPhi * cTheta);
    /*ELEM(AHRS_EKF_u, 6,0) += 0;
     ELEM(AHRS_EKF_u, 7,0) += 0;
     ELEM(AHRS_EKF_u, 8,0) += 0;
     ELEM(AHRS_EKF_u, 9,0) += 0;
     ELEM(AHRS_EKF_u, 10,0) += 0;
     ELEM(AHRS_EKF_u, 11,0) += 0;*/

    return;
}

/*-----------------------------Update with accel & gyro------------------------------------*/
void AHRS_EKF_updateAccelGyro(axis3f_t* angles, axis3f_t* velocities, axis3f_t accel, axis3f_t gyro) {
    matrix_t deltaM;
    matrixInit(&deltaM, 2, 1);
    /* Remove bias */
    float pr = gyro.x - ELEM(AHRS_EKF_u, 9, 0);
    float qr = gyro.y - ELEM(AHRS_EKF_u, 10, 0);
    float rr = gyro.z - ELEM(AHRS_EKF_u, 11, 0);

#ifdef configAHRS_EKF_CORRECT_ACCEL_OFFSET
    accel.x += ((qr * qr + rr * rr) * configAHRS_EKF_ACCEL_OFFSET_X) - (pr * qr * configAHRS_EKF_ACCEL_OFFSET_Y)
               - (pr * rr * configAHRS_EKF_ACCEL_OFFSET_Z);
    accel.y += ((pr * pr + rr * rr) * configAHRS_EKF_ACCEL_OFFSET_Y) - (pr * qr * configAHRS_EKF_ACCEL_OFFSET_X)
               - (qr * rr * configAHRS_EKF_ACCEL_OFFSET_Z);
#endif

    /* C matrix */
    //_C.zeros(); //zeros or not?
    ELEM(_C_acc, 0, 3) = -ELEM(AHRS_EKF_u, 6, 0);
    ELEM(_C_acc, 0, 4) = rr;
    ELEM(_C_acc, 0, 5) = -qr;
    ELEM(_C_acc, 0, 6) = -ELEM(AHRS_EKF_u, 3, 0);
    ELEM(_C_acc, 1, 3) = -rr;
    ELEM(_C_acc, 1, 4) = -ELEM(AHRS_EKF_u, 6, 0);
    ELEM(_C_acc, 1, 5) = pr;
    ELEM(_C_acc, 1, 6) = -ELEM(AHRS_EKF_u, 4, 0);

    /* Delta measures */
    ELEM(deltaM, 0, 0) = accel.x + ELEM(AHRS_EKF_u, 5, 0) * qr - ELEM(AHRS_EKF_u, 4, 0) * rr
                         + ELEM(AHRS_EKF_u, 6, 0) * ELEM(AHRS_EKF_u, 3, 0);
    ELEM(deltaM, 1, 0) = accel.y + ELEM(AHRS_EKF_u, 3, 0) * rr - ELEM(AHRS_EKF_u, 5, 0) * pr
                         + ELEM(AHRS_EKF_u, 6, 0) * ELEM(AHRS_EKF_u, 4, 0);

    /* Gain matrix K */
    //_M = QuadProd(_C, _P) + _R;
    QuadProd(&_C_acc, &_P, &_M);
    matrixAdd(&_M, &_R_acc, &_M);
    //_K = _P * (~_C) * (!_M);
    matrixMult_rhsT(&_P, &_C_acc, &TMP1); //TMP1 contains _P * (~_C)
    matrixInversed(&_M, &TMP2);           //TMP2 contains (!_M)
    matrixMult(&TMP1, &TMP2, &_K);

    /* Correct state vector */
    //u += _K * deltaM
    matrixMult(&_K, &deltaM, &TMP3);
    matrixAdd(&AHRS_EKF_u, &TMP3, &AHRS_EKF_u);

    /* Updated P matrix */
    //_P -= _K * _C * _P;
    //_P=(_P+(~_P))*0.5; //guarantees P to be symmetric
    matrixMult(&_K, &_C_acc, &TMP4);
    matrixMult(&TMP4, &_P, &TMP5);
    matrixSub(&_P, &TMP5, &_P);

    /* Set angles */
    angles->x = ELEM(AHRS_EKF_u, 0, 0); //u(0,0) is roll according to IMU ref. frame
    angles->y = ELEM(AHRS_EKF_u, 1, 0); //u(1,0) is pitch according to IMU ref. frame
    angles->z = ELEM(AHRS_EKF_u, 2, 0); //u(2,0) is yaw according to IMU ref. frame

    /* Set velocities */
    velocities->x = ELEM(AHRS_EKF_u, 3, 0); //u(3,0) is speed along local x axis according to IMU ref. frame
    velocities->y = ELEM(AHRS_EKF_u, 4, 0); //u(4,0) is speed along local y axis according to IMU ref. frame
    velocities->z = ELEM(AHRS_EKF_u, 5, 0); //u(5,0) is speed along local z axis according to IMU ref. frame

    matrixDelete(&deltaM);
    return;
}

//-----------------------------------Update with mag-----------------------------------------//
void AHRS_EKF_updateMag(axis3f_t* angles, axis3f_t* velocities, axis3f_t mag) {
    matrix_t deltaM, M_mag, TMP1_mag, TMP2_mag, K_mag;

    matrixInit(&deltaM, 3, 1);
    matrixInit(&M_mag, 3, 3);
    matrixInit(&TMP1_mag, 12, 3);
    matrixInit(&TMP2_mag, 3, 3);
    matrixInit(&K_mag, 12, 3);

    /* Normalize readings */
    float invNorm = 1.f / sqrtf(mag.x * mag.x + mag.y * mag.y + mag.z * mag.z);
    if (isnan(invNorm) || isinf(invNorm)) {
        invNorm = 1.f;
    }
    mag.x *= invNorm;
    mag.y *= invNorm;
    mag.z *= invNorm;

    /* Trig functions */
    float sInc = sinf(ELEM(AHRS_EKF_u, 8, 0));
    float cInc = cosf(ELEM(AHRS_EKF_u, 8, 0));
    float sPhi = sinf(ELEM(AHRS_EKF_u, 0, 0));
    float cPhi = cosf(ELEM(AHRS_EKF_u, 0, 0));
    float sTheta = sinf(ELEM(AHRS_EKF_u, 1, 0));
    float cTheta = cosf(ELEM(AHRS_EKF_u, 1, 0));
    float sPsi = sinf(ELEM(AHRS_EKF_u, 2, 0));
    float cPsi = cosf(ELEM(AHRS_EKF_u, 2, 0));

    /* C matrix */
    //_C_mag.zeros(); //zeros or not?
    float tmp1 = cPhi * sPsi - cPsi * sPhi * sTheta;
    ELEM(_C_mag, 0, 1) = -sInc * cTheta - cInc * cPsi * sTheta;
    ELEM(_C_mag, 0, 2) = -cInc * cTheta * sPsi;
    ELEM(_C_mag, 0, 7) = -cInc * sTheta - cPsi * cTheta * sInc;
    ELEM(_C_mag, 1, 0) = cInc * (sPhi * sPsi + cPhi * cPsi * sTheta) + cPhi * sInc * cTheta;
    ELEM(_C_mag, 1, 1) = cInc * cPsi * cTheta * sPhi - sInc * sPhi * sTheta;
    ELEM(_C_mag, 1, 2) = -cInc * (cPhi * cPsi + sPhi * sPsi * sTheta);
    ELEM(_C_mag, 1, 8) = sInc * tmp1 + cInc * cTheta * sPhi;
    ELEM(_C_mag, 2, 0) = cInc * tmp1 - sInc * cTheta * sPhi;
    ELEM(_C_mag, 2, 1) = cInc * cPhi * cPsi * cTheta - cPhi * sInc * sTheta;
    ELEM(_C_mag, 2, 2) = cInc * (cPsi * sPhi - cPhi * sPsi * sTheta);
    ELEM(_C_mag, 2, 8) = cInc * cPhi * cTheta - sInc * (sPhi * sPsi + cPhi * cPsi * sTheta);

    ELEM(deltaM, 0, 0) = mag.x - (cInc * cPsi * cTheta - sInc * sTheta);
    ELEM(deltaM, 1, 0) = mag.y - (sInc * cTheta * sPhi - cInc * tmp1);
    ELEM(deltaM, 2, 0) = mag.z - (cInc * (sPhi * sPsi + cPhi * cPsi * sTheta) + cPhi * sInc * cTheta);

    /* Gain matrix K */
    //_M = QuadProd(_C, _P) + _R;
    QuadProd(&_C_mag, &_P, &M_mag);
    matrixAdd(&M_mag, &_R_mag, &M_mag);
    //_K = _P * (~_C) * (!_M);
    matrixMult_rhsT(&_P, &_C_mag, &TMP1_mag); //TMP1 contains _P * (~_C)
    matrixInversed(&M_mag, &TMP2_mag);        //TMP2 contains (!_M)
    matrixMult(&TMP1_mag, &TMP2_mag, &K_mag);

    /* Correct state vector */
    //u += _K * deltaM
    matrixMult(&K_mag, &deltaM, &TMP3);
    matrixAdd(&AHRS_EKF_u, &TMP3, &AHRS_EKF_u);

    /* Updated P matrix */
    //_P -= _K * _C * _P;
    //_P=(_P+(~_P))*0.5; //guarantees P to be symmetric
    matrixMult(&K_mag, &_C_mag, &TMP4);
    matrixMult(&TMP4, &_P, &TMP5);
    matrixSub(&_P, &TMP5, &_P);

    /* Set angles */
    angles->x = ELEM(AHRS_EKF_u, 0, 0); //u(0,0) is roll according to IMU ref. frame
    angles->y = ELEM(AHRS_EKF_u, 1, 0); //u(1,0) is pitch according to IMU ref. frame
    angles->z = ELEM(AHRS_EKF_u, 2, 0); //u(2,0) is yaw according to IMU ref. frame

    /* Set velocities */
    velocities->x = ELEM(AHRS_EKF_u, 3, 0); //u(3,0) is speed along local x axis according to IMU ref. frame
    velocities->y = ELEM(AHRS_EKF_u, 4, 0); //u(4,0) is speed along local y axis according to IMU ref. frame
    velocities->z = ELEM(AHRS_EKF_u, 5, 0); //u(5,0) is speed along local z axis according to IMU ref. frame

    matrixDelete(&deltaM);
    matrixDelete(&M_mag);
    matrixDelete(&TMP1_mag);
    matrixDelete(&TMP2_mag);
    matrixDelete(&K_mag);
    return;
}

/*------------------------Update with velocity along x,y local-----------------------------*/
void AHRS_EKF_updateVelXY(axis3f_t* angles, axis3f_t* velocities, float vx, float vy, float dt_s) {
    matrix_t deltaM, R_tmp, C_tmp;

    matrixInit(&deltaM, 2, 1);
    matrixInit(&R_tmp, 2, 2);
    matrixInit(&C_tmp, 2, 12);

    /* R matrix */
    ELEM(R_tmp, 0, 0) = _r_vxy / dt_s;
    ELEM(R_tmp, 1, 1) = ELEM(R_tmp, 0, 0);

    /* C matrix */
    ELEM(C_tmp, 0, 3) = 1.f;
    ELEM(C_tmp, 1, 4) = 1.f;

    /* Delta measures */
    ELEM(deltaM, 0, 0) = vx - ELEM(AHRS_EKF_u, 3, 0);
    ELEM(deltaM, 1, 0) = vy - ELEM(AHRS_EKF_u, 4, 0);

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
    matrixAdd(&AHRS_EKF_u, &TMP3, &AHRS_EKF_u);

    /* Updated P matrix */
    //_P -= _K * C_tmp * _P;
    //_P=(_P+(~_P))*0.5; //guarantees P to be symmetric
    matrixMult(&_K, &C_tmp, &TMP4);
    matrixMult(&TMP4, &_P, &TMP5);
    matrixSub(&_P, &TMP5, &_P);

    /* Set angles */
    angles->x = ELEM(AHRS_EKF_u, 0, 0); //u(0,0) is roll according to IMU ref. frame
    angles->y = ELEM(AHRS_EKF_u, 1, 0); //u(1,0) is pitch according to IMU ref. frame
    angles->z = ELEM(AHRS_EKF_u, 2, 0); //u(2,0) is yaw according to IMU ref. frame

    /* Set velocities */
    velocities->x = ELEM(AHRS_EKF_u, 3, 0); //u(3,0) is speed along local x axis according to IMU ref. frame
    velocities->y = ELEM(AHRS_EKF_u, 4, 0); //u(4,0) is speed along local y axis according to IMU ref. frame
    velocities->z = ELEM(AHRS_EKF_u, 5, 0); //u(5,0) is speed along local z axis according to IMU ref. frame

    matrixDelete(&deltaM);
    matrixDelete(&R_tmp);
    matrixDelete(&C_tmp);

    return;
}

/*-------------------------Update with velocity along z local-------------------------------*/
void AHRS_EKF_updateVelZ(axis3f_t* angles, axis3f_t* velocities, float vz, float dt_s) {
    matrix_t C_tmp, K;

    matrixInit(&C_tmp, 1, 12);
    matrixInit(&K, 12, 1);

    /* C matrix */
    ELEM(C_tmp, 0, 5) = 1.f;

    /* Delta measures */
    float deltaM = vz - ELEM(AHRS_EKF_u, 5, 0);

    /* Gain matrix K (change K with _K in correction and update of P) */
    /*_M = QuadProd(C_tmp,_P) + (_r_vz / dt_s);
   _K = _P * (~C_tmp) * (!_M);*/

    /* Faster Gain matrix K */
    float inv_m = 1.f / (ELEM(_P, 5, 5) + (_r_vz / dt_s));
    ELEM(K, 0, 0) = ELEM(_P, 0, 5) * inv_m;
    ELEM(K, 1, 0) = ELEM(_P, 1, 5) * inv_m;
    ELEM(K, 2, 0) = ELEM(_P, 2, 5) * inv_m;
    ELEM(K, 3, 0) = ELEM(_P, 3, 5) * inv_m;
    ELEM(K, 4, 0) = ELEM(_P, 4, 5) * inv_m;
    ELEM(K, 5, 0) = ELEM(_P, 5, 5) * inv_m;
    ELEM(K, 6, 0) = ELEM(_P, 6, 5) * inv_m;
    ELEM(K, 7, 0) = ELEM(_P, 7, 5) * inv_m;
    ELEM(K, 8, 0) = ELEM(_P, 8, 5) * inv_m;
    ELEM(K, 9, 0) = ELEM(_P, 9, 5) * inv_m;
    ELEM(K, 10, 0) = ELEM(_P, 10, 5) * inv_m;
    ELEM(K, 11, 0) = ELEM(_P, 11, 5) * inv_m;

    /* Correct state vector */
    //u += K * deltaM;
    matrixMultScalar(&K, deltaM, &TMP3);
    matrixAdd(&AHRS_EKF_u, &TMP3, &AHRS_EKF_u);

    /* Updated P matrix */
    //_P -= K * C_tmp * _P;
    //_P=(_P+(~_P))*0.5; //guarantees P to be symmetric
    matrixMult(&K, &C_tmp, &TMP4);
    matrixMult(&TMP4, &_P, &TMP5);
    matrixSub(&_P, &TMP5, &_P);

    /* Set angles */
    angles->x = ELEM(AHRS_EKF_u, 0, 0); //u(0,0) is roll according to IMU ref. frame
    angles->y = ELEM(AHRS_EKF_u, 1, 0); //u(1,0) is pitch according to IMU ref. frame
    angles->z = ELEM(AHRS_EKF_u, 2, 0); //u(2,0) is yaw according to IMU ref. frame

    /* Set velocities */
    velocities->x = ELEM(AHRS_EKF_u, 3, 0); //u(3,0) is speed along local x axis according to IMU ref. frame
    velocities->y = ELEM(AHRS_EKF_u, 4, 0); //u(4,0) is speed along local y axis according to IMU ref. frame
    velocities->z = ELEM(AHRS_EKF_u, 5, 0); //u(5,0) is speed along local z axis according to IMU ref. frame

    matrixDelete(&C_tmp);
    matrixDelete(&K);

    return;
}

//-------------------------Update with velocity along N,E global------------------------------//
void AHRS_EKF_updateVelNE(axis3f_t* angles, axis3f_t* velocities, float vN, float vE, float dt_s) {
    matrix_t deltaM, R_tmp, C_tmp;

    matrixInit(&deltaM, 2, 1);
    matrixInit(&R_tmp, 2, 2);
    matrixInit(&C_tmp, 2, 12);

    /* R matrix */
    ELEM(R_tmp, 0, 0) = _r_vne / dt_s;
    ELEM(R_tmp, 1, 1) = ELEM(R_tmp, 0, 0);

    /* Trig functions */
    float sPhi = sinf(ELEM(AHRS_EKF_u, 0, 0));
    float cPhi = cosf(ELEM(AHRS_EKF_u, 0, 0));
    float sTheta = sinf(ELEM(AHRS_EKF_u, 1, 0));
    float cTheta = cosf(ELEM(AHRS_EKF_u, 1, 0));
    float sPsi = sinf(ELEM(AHRS_EKF_u, 2, 0));
    float cPsi = cosf(ELEM(AHRS_EKF_u, 2, 0));

    /* C matrix */
    //C_tmp.zeros(); //zeros or not?
    float tmp1 = sPhi * sPsi + cPhi * cPsi * sTheta;
    float tmp2 = cPhi * cPsi + sPhi * sPsi * sTheta;
    float tmp3 = cPhi * sPsi - cPsi * sPhi * sTheta;
    float tmp4 = cPsi * sPhi - cPhi * sPsi * sTheta;
    ELEM(C_tmp, 0, 0) = ELEM(AHRS_EKF_u, 4, 0) * tmp1 + ELEM(AHRS_EKF_u, 5, 0) * tmp3;
    ELEM(C_tmp, 0, 1) = ELEM(AHRS_EKF_u, 5, 0) * cPhi * cPsi * cTheta - ELEM(AHRS_EKF_u, 3, 0) * cPsi * sTheta
                        + ELEM(AHRS_EKF_u, 4, 0) * cPsi * cTheta * sPhi;
    ELEM(C_tmp, 0, 2) =
        ELEM(AHRS_EKF_u, 5, 0) * -ELEM(AHRS_EKF_u, 4, 0) * tmp2 - ELEM(AHRS_EKF_u, 3, 0) * cTheta * sPsi;
    ELEM(C_tmp, 0, 3) = cPsi * cTheta;
    ELEM(C_tmp, 0, 4) = -tmp3;
    ELEM(C_tmp, 0, 5) = tmp1;
    ELEM(C_tmp, 1, 0) = -ELEM(AHRS_EKF_u, 4, 0) * tmp4 - ELEM(AHRS_EKF_u, 5, 0) * tmp2;
    ELEM(C_tmp, 1, 1) = ELEM(AHRS_EKF_u, 5, 0) * cPhi * cTheta * sPsi - ELEM(AHRS_EKF_u, 3, 0) * sPsi * sTheta
                        + ELEM(AHRS_EKF_u, 4, 0) * cTheta * sPhi * sPsi;
    ELEM(C_tmp, 1, 2) =
        ELEM(AHRS_EKF_u, 5, 0) * tmp1 - ELEM(AHRS_EKF_u, 4, 0) * tmp3 + ELEM(AHRS_EKF_u, 3, 0) * cPsi * cTheta;
    ELEM(C_tmp, 1, 3) = cTheta * sPsi;
    ELEM(C_tmp, 1, 4) = tmp2;
    ELEM(C_tmp, 1, 5) = tmp4;

    /* Delta measures */
    ELEM(deltaM, 0, 0) =
        vN - (ELEM(AHRS_EKF_u, 5, 0) * tmp1 - ELEM(AHRS_EKF_u, 4, 0) * tmp3 + ELEM(AHRS_EKF_u, 3, 0) * cPsi * cTheta);
    ELEM(deltaM, 1, 0) =
        vE - (ELEM(AHRS_EKF_u, 4, 0) * tmp2 - ELEM(AHRS_EKF_u, 5, 0) * tmp4 + ELEM(AHRS_EKF_u, 3, 0) * cTheta * sPsi);

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
    matrixAdd(&AHRS_EKF_u, &TMP3, &AHRS_EKF_u);

    /* Updated P matrix */
    //_P -= _K * C_tmp * _P;
    //_P=(_P+(~_P))*0.5; //guarantees P to be symmetric
    matrixMult(&_K, &C_tmp, &TMP4);
    matrixMult(&TMP4, &_P, &TMP5);
    matrixSub(&_P, &TMP5, &_P);

    /* Set angles */
    angles->x = ELEM(AHRS_EKF_u, 0, 0); //u(0,0) is roll according to IMU ref. frame
    angles->y = ELEM(AHRS_EKF_u, 1, 0); //u(1,0) is pitch according to IMU ref. frame
    angles->z = ELEM(AHRS_EKF_u, 2, 0); //u(2,0) is yaw according to IMU ref. frame

    /* Set velocities */
    velocities->x = ELEM(AHRS_EKF_u, 3, 0); //u(3,0) is speed along local x axis according to IMU ref. frame
    velocities->y = ELEM(AHRS_EKF_u, 4, 0); //u(4,0) is speed along local y axis according to IMU ref. frame
    velocities->z = ELEM(AHRS_EKF_u, 5, 0); //u(5,0) is speed along local z axis according to IMU ref. frame

    matrixDelete(&deltaM);
    matrixDelete(&R_tmp);
    matrixDelete(&C_tmp);

    return;
}

//--------------------------Update with velocity along D global--------------------------------//
void AHRS_EKF_updateVelD(axis3f_t* angles, axis3f_t* velocities, float vD, float dt_s) {
    matrix_t C_tmp, K, M;

    matrixInit(&C_tmp, 1, 12);
    matrixInit(&K, 12, 1);
    matrixInit(&M, 1, 1);

    /* Trig functions */
    float sPhi = sinf(ELEM(AHRS_EKF_u, 0, 0));
    float cPhi = cosf(ELEM(AHRS_EKF_u, 0, 0));
    float sTheta = sinf(ELEM(AHRS_EKF_u, 1, 0));
    float cTheta = cosf(ELEM(AHRS_EKF_u, 1, 0));

    /* C matrix */
    //C_tmp.zeros(); //zeros or not?
    ELEM(C_tmp, 0, 0) = ELEM(AHRS_EKF_u, 4, 0) * cPhi * cTheta - ELEM(AHRS_EKF_u, 5, 0) * cTheta * sPhi;
    ELEM(C_tmp, 0, 1) = -ELEM(AHRS_EKF_u, 3, 0) * cTheta - ELEM(AHRS_EKF_u, 5, 0) * cPhi * sTheta
                        - ELEM(AHRS_EKF_u, 4, 0) * sPhi * sTheta;
    ELEM(C_tmp, 0, 3) = -sTheta;
    ELEM(C_tmp, 0, 4) = cTheta * sPhi;
    ELEM(C_tmp, 0, 5) = cPhi * cTheta;

    /* Delta measures */
    float deltaM = vD
                   - (ELEM(AHRS_EKF_u, 5, 0) * cPhi * cTheta - ELEM(AHRS_EKF_u, 3, 0) * sTheta
                      + ELEM(AHRS_EKF_u, 4, 0) * cTheta * sPhi);

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
    matrixAdd(&AHRS_EKF_u, &TMP3, &AHRS_EKF_u);

    /* Updated P matrix */
    //_P -= _K * C_tmp * _P;
    //_P=(_P+(~_P))*0.5; //guarantees P to be symmetric
    matrixMult(&K, &C_tmp, &TMP4);
    matrixMult(&TMP4, &_P, &TMP5);
    matrixSub(&_P, &TMP5, &_P);

    /* Set angles */
    angles->x = ELEM(AHRS_EKF_u, 0, 0); //u(0,0) is roll according to IMU ref. frame
    angles->y = ELEM(AHRS_EKF_u, 1, 0); //u(1,0) is pitch according to IMU ref. frame
    angles->z = ELEM(AHRS_EKF_u, 2, 0); //u(2,0) is yaw according to IMU ref. frame

    /* Set velocities */
    velocities->x = ELEM(AHRS_EKF_u, 3, 0); //u(3,0) is speed along local x axis according to IMU ref. frame
    velocities->y = ELEM(AHRS_EKF_u, 4, 0); //u(4,0) is speed along local y axis according to IMU ref. frame
    velocities->z = ELEM(AHRS_EKF_u, 5, 0); //u(5,0) is speed along local z axis according to IMU ref. frame

    return;
}

//-----------------------------------Starting values-----------------------------------------//
void AHRS_EKF_reset(axis3f_t* angles, axis3f_t* velocities, float phi0, float theta0, float psi0) {
    matrixSet(&AHRS_EKF_u, 0, 0, phi0);
    matrixSet(&AHRS_EKF_u, 1, 0, theta0);
    matrixSet(&AHRS_EKF_u, 2, 0, psi0);
    matrixSet(&AHRS_EKF_u, 3, 0, 0);
    matrixSet(&AHRS_EKF_u, 4, 0, 0);
    matrixSet(&AHRS_EKF_u, 5, 0, 0);
    // matrixSet(&AHRS_EKF_u, 6, 0, c_damp_0);
    // matrixSet(&AHRS_EKF_u, 8, 0, incl_0);

    /* Set angles */
    angles->x = ELEM(AHRS_EKF_u, 0, 0); //u(0,0) is roll according to IMU ref. frame
    angles->y = ELEM(AHRS_EKF_u, 1, 0); //u(1,0) is pitch according to IMU ref. frame
    angles->z = ELEM(AHRS_EKF_u, 2, 0); //u(2,0) is yaw according to IMU ref. frame

    /* Set velocities */
    velocities->x = ELEM(AHRS_EKF_u, 3, 0); //u(3,0) is speed along local x axis according to IMU ref. frame
    velocities->y = ELEM(AHRS_EKF_u, 4, 0); //u(4,0) is speed along local y axis according to IMU ref. frame
    velocities->z = ELEM(AHRS_EKF_u, 5, 0); //u(5,0) is speed along local z axis according to IMU ref. frame

    return;
}

//-------------------------------------Input noises------------------------------------------//
void AHRS_EKF_setInputNoises(float gxy, float gz, float az, float c_damp, float b_az, float incl, float b_g) {
    matrixSet(&_W, 0, 0, gxy * configAHRS_EKF_LOOP_TIME_S);
    matrixSet(&_W, 1, 1, gxy * configAHRS_EKF_LOOP_TIME_S);
    matrixSet(&_W, 2, 2, gz * configAHRS_EKF_LOOP_TIME_S);
    matrixSet(&_W, 3, 3, az * configAHRS_EKF_LOOP_TIME_S);
    matrixSet(&_W, 4, 4, c_damp * configAHRS_EKF_LOOP_TIME_S);
    matrixSet(&_W, 5, 5, b_az * configAHRS_EKF_LOOP_TIME_S);
    matrixSet(&_W, 6, 6, incl * configAHRS_EKF_LOOP_TIME_S);
    matrixSet(&_W, 7, 7, b_g * configAHRS_EKF_LOOP_TIME_S);
    matrixSet(&_W, 8, 8, b_g * configAHRS_EKF_LOOP_TIME_S);
    matrixSet(&_W, 9, 9, b_g * configAHRS_EKF_LOOP_TIME_S);

    return;
}

//----------------------------------Output accel noises--------------------------------------//
void AHRS_EKF_setAccelNoise(float axy) {
    float inv_loop_time_pred_update_acc_s = 1.0f / configAHRS_EKF_LOOP_TIME_S;
    matrixSet(&_R_acc, 0, 0, axy * inv_loop_time_pred_update_acc_s);
    matrixSet(&_R_acc, 1, 1, axy * inv_loop_time_pred_update_acc_s);
    return;
}

//-----------------------------------Output mag noises---------------------------------------//
void AHRS_EKF_setMagNoise(float m) {
    float inv_loop_time_mag_s = 1.0f / configAHRS_EKF_MAG_LOOP_TIME_S;
    matrixSet(&_R_mag, 0, 0, m * inv_loop_time_mag_s);
    matrixSet(&_R_mag, 1, 1, m * inv_loop_time_mag_s);
    matrixSet(&_R_mag, 2, 2, m * inv_loop_time_mag_s);
    return;
}

/*----------------------------Get values from state vector----------------------------------*/
float AHRS_EKF_getStateValue(uint8_t idx) { return matrixGet(&AHRS_EKF_u, idx, 0); }
