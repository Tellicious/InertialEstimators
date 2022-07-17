//
//  main.cpp
//  Test_Kalman
//
//  Created by Andrea Vivani on 22/4/15.
//  Copyright (c) 2015 Andrea Vivani. All rights reserved.
//

//
// CF ha KP=1 e Ki=0.05

#include <iostream>
#include <stdio.h>
#include <chrono>
#include "MatLib.h"
#include "AHRS_EKF_AV.h"
#include "IMU_EKF_AV.h"
#include "AHRS_Madgwick.h"
#include "IMU_Madgwick.h"
#include "AHRS_Attitude_EKF.h"
#include "AHRS_Attitude_SO3.h"
#include "Utilities.h"
#include "attitude_estimator_so3_main.h"

extern "C"{
#include "AttitudeEKF.h"
#include "IMU_EKF.h"
}

template <typename T> void mprint(MatrixX<T>& r) {
    std::cout << "\n" << (int) r.rows() << "x" << (int) r.columns() << "\n";
    for(int i=0; i<r.rows(); i++) {
        for (int j=0; j<r.columns(); j++)
            printf("%4.6f\t",(double) r(i,j));
        std::cout << "\n";
    }
};

template <typename T> void mprint_file(FILE *myfile, MatrixX<T>& r) {
    for(int i=0; i<r.rows(); i++) {
        for (int j=0; j<r.columns(); j++)
            fprintf(myfile,"%4.6f,",(double) r(i,j));
        fprintf(myfile, "\n");
    }
};

void read_csv_all(FILE *myfile, uint16_t rows, MatrixXf &t_s, MatrixXf &g, MatrixXf &a, MatrixXf &m, MatrixXf &v, MatrixXf &V){
    for(uint16_t ii=0;ii<rows;ii++){
        fscanf(myfile,"%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", &t_s(ii,0), &g(ii,0), &g(ii,1), &g(ii,2), &a(ii,0), &a(ii,1), &a(ii,2), &m(ii,0), &m(ii,1), &m(ii,2), &v(ii,0), &v(ii,1), &v(ii,2), &V(ii,0), &V(ii,1), &V(ii,2));
    }
    return;
};

void read_csv_short(FILE *myfile, uint16_t rows, MatrixXf &t_s, MatrixXf &m){
    for(uint16_t ii=0;ii<rows;ii++){
        fscanf(myfile,"%f,%f,%f,%f\n", &t_s(ii,0), &m(ii,0), &m(ii,1), &m(ii,2));
    }
    return;
};

// CLASSES DEFINITION

AHRS_EKF_AV AV_AHRS(9.81, 5e-3);
IMU_EKF_AV AV_IMU(9.81, 5e-3);
AHRS_Madgwick M_AHRS;
IMU_Madgwick M_IMU;
AHRS_Attitude_EKF PX4_AHRS(9.81, 5e-3);
AHRS_Attitude_SO3 PX4_CF(5e-3);

AHRS_State_t ahrs_tmp;

int main(int argc, const char * argv[]) {
    AttitudeEKF_initialize();
    PX4_AHRS.set_inclination();
    IMU_EKF_init(&ahrs_tmp, 9.81, 5e-3);
    
    // READ DATA FROM CSV FILES
    FILE * data1;
    FILE * data2;
    MatrixXf t_samp(4001,1);
    MatrixXf gyro(4001,3);
    MatrixXf acc(4001,3);
    MatrixXf mag1(4001,3);
    MatrixXf v(4001,3);
    MatrixXf V(4001,3);
    MatrixXf t_samp2(501,1);
    MatrixXf mag2(501,3);
    data1=fopen("long.csv", "r");
    data2=fopen("short.csv", "r");
    read_csv_all(data1, 4001, t_samp, gyro, acc, mag1, v, V);
    read_csv_short(data2, 501, t_samp2, mag2);
    
    // OUTPUT FILES
    FILE * VivaniAHRS;
    FILE * VivaniIMU;
    FILE * VivaniIMU_C;
    FILE * MadgwickAHRS;
    FILE * MadgwickIMU;
    FILE * PX4AHRS;
    FILE * PX4AHRSC;
    FILE * PX4CF;
    
    VivaniAHRS=fopen("VivaniAHRS.csv", "w");
    MatrixXf printU=~AV_AHRS.u;
    mprint_file(VivaniAHRS, printU);
    
    VivaniIMU=fopen("VivaniIMU.csv", "w");
    printU=~AV_IMU.u;
    mprint_file(VivaniIMU, printU);
    
    VivaniIMU_C=fopen("VivaniIMU_C.csv", "w");
    fprintf(VivaniIMU_C, "%4.6f,%4.6f,%4.6f,%4.6f,%4.6f,%4.6f,%4.6f\n",IMU_EKF_u.data[0],IMU_EKF_u.data[1],IMU_EKF_u.data[2],IMU_EKF_u.data[3],IMU_EKF_u.data[4],IMU_EKF_u.data[5],IMU_EKF_u.data[6]);
    
    MadgwickAHRS=fopen("MadgwickAHRS.csv","w");
    fprintf(MadgwickAHRS, "%4.6f,%4.6f,%4.6f\n",M_AHRS.Roll,M_AHRS.Pitch,M_AHRS.Yaw);
    
    MadgwickIMU=fopen("MadgwickIMU.csv","w");
    fprintf(MadgwickIMU, "%4.6f,%4.6f\n",M_IMU.Roll,M_IMU.Pitch);
    
    PX4AHRS=fopen("PX4AHRS.csv","w");
    fprintf(PX4AHRS, "0,0,0\n");
    
    PX4AHRSC=fopen("PX4_C.csv","w");
    fprintf(PX4AHRSC,"%4.6f,%4.6f,%4.6f\n",0.f, 0.f, 0.f);
    
    PX4CF=fopen("PX4_CF.csv","w");
    fprintf(PX4CF,"%4.6f,%4.6f,%4.6f\n",0.f, 0.f, 0.f);
    
    // INITIALIZATION
    float gx_old=0;
    float gy_old=0;
    float gz_old=0;
    float az_old=0;
    uint16_t ind_mag=0;
    
    // LOOP
    for (uint16_t ii=0; ii<t_samp.rows()-1; ii++) {
        
        // AHRS Andrea Vivani
        AV_AHRS.prediction(gx_old, gy_old, gz_old, az_old);
        AV_AHRS.update_accel(gyro(ii,0), gyro(ii,1), gyro(ii,2), acc(ii,0), acc(ii,1));
        AV_AHRS.update_mag(mag1(ii,0), mag1(ii,1), mag1(ii,2));
        //AV_AHRS.update_vel_xy(v(ii,0), v(ii,1), 5e-3);
        //AV_AHRS.update_vel_z(v(ii,2), 5e-3);
        AV_AHRS.update_vel_ne(V(ii,0), V(ii,1), 5e-3);
        AV_AHRS.update_vel_d(V(ii,2), 5e-3);
        /*if (t_samp(ii,0)>=t_samp2(ind_mag,0)){
         AV_AHRS.update_mag(mag2(ind_mag,0), mag2(ind_mag,1), mag2(ind_mag,2));
         ind_mag++;
         }*/
        
        
        // IMU Andrea Vivani
        AV_IMU.prediction(gx_old, gy_old, gz_old, az_old);
        AV_IMU.update(gyro(ii,0), gyro(ii,1), gyro(ii,2), acc(ii,0), acc(ii,1));
        //AV_IMU.update_vel_xy(v(ii,0), v(ii,1), 5e-3);
        //AV_IMU.update_vel_z(v(ii,2), 5e-3);
        //AV_IMU.update_vel_d(V(ii,2), 5e-3);
        
        IMU_EKF_prediction(gx_old, gy_old, gz_old, az_old);
        IMU_EKF_update(&ahrs_tmp, gyro(ii,0), gyro(ii,1), gyro(ii,2), acc(ii,0), acc(ii,1));
        //IMU_EKF_update_vel_xy(&ahrs_tmp, v(ii,0), v(ii,1), 5e-3);
        //IMU_EKF_update_vel_z(&ahrs_tmp, v(ii,2), 5e-3); //this one is different!!!
        //IMU_EKF_update_vel_d(&ahrs_tmp, V(ii,2), 5e-3); //this one is different
        
        // AHRS Sebastian Madgwick
        M_AHRS.compute(gyro(ii,0), gyro(ii,1), gyro(ii,2), acc(ii,0), acc(ii,1), acc(ii,2), mag1(ii,0), mag1(ii,1), mag1(ii,2), 5e-3);
        
        // IMU Sebastian Madgwick
        M_IMU.compute(gyro(ii,0), gyro(ii,1), gyro(ii,2), acc(ii,0), acc(ii,1), acc(ii,2), 5e-3);
        
        // AHRS PIXHAWK (EKF)
        PX4_AHRS.prediction();
        PX4_AHRS.update_gyro(gyro(ii,0), gyro(ii,1), gyro(ii,2));
        PX4_AHRS.update_accel(acc(ii,0), acc(ii,1), acc(ii,2));
        PX4_AHRS.update_mag(mag1(ii,0), mag1(ii,1), mag1(ii,2));
        PX4_AHRS.euler_angles();
        
        // ORIGINAL AHRS PIXHAWK (EKF, C CODE)
        unsigned char m_disp[3]={1,1,1};
        float z[]={gyro(ii,0), gyro(ii,1), gyro(ii,2), acc(ii,0), acc(ii,1), acc(ii,2), mag1(ii,0), mag1(ii,1), mag1(ii,2) };
        float J[]={1, 0, 0, 0, 1, 0, 0, 0, 1};
        float xa_apo[12];
        float Pa_apo[144];
        float Rot_matrix[9];
        float Angles[3];
        float debug[4];
        AttitudeEKF(1, 0, m_disp, 5e-3, z, 1e-4, 0.08, 0.009, 0.005, 0.0008, 10000, 100, J, xa_apo, Pa_apo, Rot_matrix, Angles, debug);
        
        // AHRS PIXHAWK COMPLEMENTARY FILTER
        PX4_CF.compute(gyro(ii,0), gyro(ii,1), gyro(ii,2), acc(ii,0), acc(ii,1), acc(ii,2), mag1(ii,0), mag1(ii,1), mag1(ii,2));
        
        NonlinearSO3AHRSupdate(gyro(ii,0), gyro(ii,1), gyro(ii,2), acc(ii,0), acc(ii,1), acc(ii,2), mag1(ii,0), mag1(ii,1), mag1(ii,2), 1, 0.05, 5e-3);
        
        // TEMPORARY VARIABLES
        gx_old=gyro(ii,0);
        gy_old=gyro(ii,1);
        gz_old=gyro(ii,2);
        az_old=acc(ii,2);
        
        // PRINT TO FILE
        MatrixXf printU=~AV_AHRS.u;
        mprint_file(VivaniAHRS, printU);
        printU=~AV_IMU.u;
        mprint_file(VivaniIMU, printU);
        fprintf(VivaniIMU_C, "%4.6f,%4.6f,%4.6f,%4.6f,%4.6f,%4.6f,%4.6f\n",IMU_EKF_u.data[0],IMU_EKF_u.data[1],IMU_EKF_u.data[2],IMU_EKF_u.data[3],IMU_EKF_u.data[4],IMU_EKF_u.data[5],IMU_EKF_u.data[6]);
        fprintf(MadgwickAHRS, "%4.6f,%4.6f,%4.6f\n",M_AHRS.Roll,M_AHRS.Pitch,M_AHRS.Yaw);
        fprintf(MadgwickIMU, "%4.6f,%4.6f\n",M_IMU.Roll,M_IMU.Pitch);
        fprintf(PX4AHRS,"%4.6f,%4.6f,%4.6f\n",PX4_AHRS.Roll, PX4_AHRS.Pitch, PX4_AHRS.Yaw);
        fprintf(PX4AHRSC,"%4.6f,%4.6f,%4.6f\n",Angles[0], Angles[1], Angles[2]);
        //fprintf(PX4CF,"%4.6f,%4.6f,%4.6f\n",PX4_CF.Roll, PX4_CF.Pitch, PX4_CF.Yaw);
        fprintf(PX4CF,"%4.6f,%4.6f,%4.6f\n",euler[0],euler[1],euler[2]);
    }
    
    // CLOSE ALL PREVIOUSLY OPENED FILES
    fclose(VivaniAHRS);
    fclose(VivaniIMU);
    fclose(VivaniIMU_C);
    fclose(MadgwickAHRS);
    fclose(MadgwickIMU);
    fclose(PX4AHRS);
    fclose(PX4AHRSC);
    fclose(PX4CF);
    return 0;
}
