# Changelog

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