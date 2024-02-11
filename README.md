<p align="center"> <img src="https://github.com/Tellicious/ADV-utils/assets/9076397/3ec512f1-2de6-4226-bc07-e4bfdd686a28" width=50% height=50%> </p>

# Collection of intertial state estimators (attitude + heading)

## C libraries included:
### IMU
#### Attitude-only estimation based on accelerometer and gyroscope 
- ***IMU_EKF:*** quadcopter-specific EKF. Features accelerometer CoG offset compensation and includes optional correction via vertical speed reading (e.g. from `altitudeKF`)
- ***IMU_quaternionST:*** quaternion-based filter based on STEVAL-FCU001v1 code

### Altitude and Rate-of-Climb
- ***altitudeKF:*** Kalman filter based on accelerometer, barometer and ToF (optional) readings

## C++ libraries included:
### AHRS
#### Attitude and heading estimation based on accelerometer, gyroscope and magnetometer
- ***AHRS_Attitude_EKF:*** EKF based on PX4 code
- ***AHRS_EKF_AV:*** quadcopter-specific EKF
- ***AHRS_Madgwick:*** quaternion-based filter based on Sebastian Madgwick work
- ***AHRS_Attitude_SO3:*** quaternion-based complementary filter based on Robert Mahony work

### IMU
#### Attitude-only estimation based on accelerometer and gyroscope 
- ***IMU_EKF_AV:*** quadcopter-specific EKF (same algorithm as `IMU_EKF` C library)
- ***IMU_Madgwick:*** quaternion-based filter based on Sebastian Madgwick work

## Warning
*C and C++ libraries, despite being based on the same algorithms, can contain different features*

*All C libraries require `ADV-utils` to run*

*All C++ libraries require `MatLib` to run*
