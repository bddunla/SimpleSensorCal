/*
 * File: magCal.h
 *
 * 
 * Description:
 * 
 */

#ifndef __MAG_CAL_H__
#define __MAG_CAL_H__

#include <stdint>

#ifdef __cplusplus
extern "C" {
#endif

// Calibration data. Applied to sensor data as follows:
// mag_calibrated_sample = scaleRotate * (mag_raw - offset)
typedef struct {
    float offsetX;
    float offsetY;
    float[4] scaleRotate; // 2x2 matrix in row-major order]
} mag_calibration_t;

void magCalDualBuffEllipseFit(float *X_data, float *Y_data, uint32_t size);

#ifdef __cplusplus
}
#endif

#endif // __MAG_CAL_H__