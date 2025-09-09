
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <vector>
#include "magCal.h"
#include <cmath>

#define DATA_NUMS 100

/* Prototypes */
void magCalEigenEllipseFit(Eigen::MatrixXf X_data, Eigen::MatrixXf Y_data, uint32_t size);

/*
 * Function: magCalEllipseFit
 *
 * Description: Takes in two arrays of floats, X and Y, for fitting
 *              sensor characterization elipse to.
 * Note: Input data is highly sensitive to clustering. For best
 *       results, input data evenly covers full "circle" of
 *       sensor's collection domain. Either from collection
 *       technique or preprocessing
 * Params: X_data - X coordinates from each sample, alligned with Y_data
 *         Y_data - Y coordinates from each sample, alligned with X_data
 *         size   - number of (X,Y) samples, i.e. length of input arrays
 */
void magCalDualBuffEllipseFit(mag_calibration_t *calData, float *X_data, float *Y_data, uint32_t size)
{
    Eigen::Map<Eigen::MatrixXf> X_mat(X_data, size, 1);
    Eigen::Map<Eigen::MatrixXf> Y_mat(Y_data, size, 1);
    magCalEigenEllipseFit(calData, X_mat, Y_mat, size);
}


void magCalEigenEllipseFit(mag_calibration_t *calData, Eigen::MatrixXf X_data, Eigen::MatrixXf Y_data, uint32_t size)
{
    // Form A matrix of 2 dimentional input variables, where Ax = b when, x is the general ellipse coefficients
    Eigen::MatrixXf matA;
    matA.resize(size,5);
    matA.col(0) = X_data.array()*Y_data.array();
    matA.col(1) = Y_data.array().square();
    matA.col(2) = X_data.array();
    matA.col(3) = Y_data.array();
    matA.col(4) = Eigen::MatrixXf::Ones(matA.rows(),1).array();

    // Solve for the least squares solution to vector x.
    Eigen::MatrixXf vecb = X_data.cwiseProduct(X_data);
    Eigen::MatrixXf x_hat = (matA.transpose() * matA).ldlt().solve(matA.transpose() * vecb);

    // Separate out the coefficients
    float A = -1.0f;
    float B = x_hat(0);
    float C = x_hat(1);
    float D = x_hat(2);
    float E = x_hat(3);
    float F = x_hat(4);

    /*
     * Next the general ellipse found is converted into 6 calibration values.
     * 2 numbers represent the translation offset of the ellipse center
     * 4 represent a 2x2 matrix to scale and rotate.
     * 
     * Although translation, rotate, and scale can be represented with a 3x3 matrix
     * this will likely be slower on MCU processors.
     */

    // Ellipse center offset
    float denominator = B*B - 4*A*C;
    float offset[2];
    offset[0] = (2*C*D - B*E) / denominator;
    offset[1] = (2*A*E - B*D) / denominator;

    // Rotate and scale
    float numeratorLeft = 2 * A*E*E + C*D*D - B*D*E + denominator*F;
    float numeratorRightX = (A + C) + std::sqrt((A - C) * (A - C) + B * B);
    float numeratorRightY = (A + C) - std::sqrt((A - C) * (A - C) + B * B);
    float numeratorX = numeratorLeft * numeratorRightX;
    float numeratorY = numeratorLeft * numeratorRightY;

    xScale = -std::sqrt(numeratorX)/denominator;
    yScale = -std::sqrt(numeratorY)/denominator;

    // Compute rotation angle
    angle = 0.5f * std::atan2(-B, C - A);

    // Build scale and rotate matrix
    Eigen::Matrix2f scaleRotate;

    scaleRotate << scaleX * std::cos(theta), -scaleY * std::sin(theta),
                   scaleX * std::sin(theta),  scaleY * std::cos(theta);


    // Apply rotate and scale then offset in forward direction.
    // Thus, reversing the sensor data offsets, then multiplies by inverse scale and rotate
    
    scaleRotate = scaleRotate.inverse();

    // Output calibration values
    calData->offsetX = offset[0];
    calData->offsetY = offset[1];
    calData->scaleRotate[0] = scaleRotate(0,0);
    calData->scaleRotate[1] = scaleRotate(0,1);
    calData->scaleRotate[2] = scaleRotate(1,0);
    calData->scaleRotate[3] = scaleRotate(1,1);
}