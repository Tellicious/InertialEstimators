//  InertialEstimator.h
//
//
//  Created by Andrea Vivani on 15/4/15.
//  Copyright (c) 2015 Andrea Vivani. All rights reserved.
//
// AXES DIRECTIONS (SENSORS AND ESTIMATION): X POINTING FORWARD (ROLL, PHI), Y POINTING RIGHT (PITCH, THETA), Z POINTING DOWN (YAW, PSI)
// IF SENSORS HAVE DIFFERENT ORIENTATION, ROTATE READINGS PRIOR THAN INPUT THEM TO THE ESTIMATORS

#ifndef InertialEstimator_H_
#define InertialEstimator_H_

#include "IMU_Madgwick.h"
#include "AHRS_Madgwick.h"
#include "IMU_EKF_AV.h"
#include "AHRS_EKF_AV.h"
#include "AHRS_Attitude_EKF.h"
#include "AHRS_Attitude_SO3.h"

#endif