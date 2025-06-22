# Changelog

## v1.3.2

**Bugfix:**
- Formatting fix in `altitudeKF`

## v1.3.1

**Bugfix:**
- Corrected a bug in `altitudeKF` LIDAR correction by adding a missing `fasbf`

## v1.3.0

**New features:**
- Added Kalman filter gain calculation during initialization to `altitudeKF`, using the `DARE` function from ADV-Utils v1.13.0
- Added ground effect detection and correction to `altitudeKF` by defining `configALTITUDE_KF_DETECT_GROUND_EFFECT`

## v1.2.3

**Bugfix:**
- Added a missing _#endif_ in `altitudeKF`
- Fixed naming of filter gains

## v1.2.2

**Bugfix:**
- Corrected `altitudeKF` MATLAB code to show gain L instead of M
- Corrected `altitudeKF` functions to always use predicted value when multiple corrections are occurring

## v1.2.1

**Improvements:**
- Removed faster math defines as already included in `basicMath` from ADVUtils
- Changed filter initialization in `altitudeKF` to use the new `IIRFilterInitHP` method

**Bugfix:**
- Removed remaining calls to `fastInvSqrt`

## v1.2.0

**Improvements:**
- Added possibility of using approximate altitude calculation based on ISA conditions (t_SL = 15°C, QNH = 101325 Pa) to `altitudeKF` by defining `configALTITUDE_KF_USE_APPROX_ALTITUDE` 

## v1.1.1

**Improvements:**
- Added possibility of using faster math functions also to `IMU_quaternionST`
- Added accelerometer offset correction also to `IMU_quaternionST`
- Other minor improvements to `IMU_quaternionST`

## v1.1.0

**New Features:**
- Added possibility of using faster math functions by adding `USE_FAST_MATH` to compile definitions, swapping all `sinf`, `cosf` and `sqrtf` occurrences with faster variants `fastSin`, `fastCos`, `fastSqrt` and `fastInvSqrt` from `basicMath`

## v1.0.6

**Bugfix:**
- Typo correction in `altitudeKF.c` that was stopping compilation on case-sensitive filesystems

## v1.0.5

**Bugfix:**
- Fixed yaw calculation logic in `AHRS_EKF.c` so that now it always returns a value from 0 to 2 * pi
- Fixed sensors and axes orientation in`AHRS_Madgwick` and `IMU_Madgwick`
- Fixed GitHub action

## v1.0.4

**Bugfix:**
- Changed yaw calculation logic in `AHRS_EKF.c` so that now it is always returned from 0 to pi

## v1.0.3

**Bugfix:**
- Fixed a bug in `AHRS_EKF.c`

## v1.0.2

**Improvements:**
- Removed useless comment in `IMU_EKF.h`

## v1.0.1

**Improvements:**
- Minor improvement in `IMU_Madgwick` and `AHRS_Madgwick`

## v1.0.0

- Initial release