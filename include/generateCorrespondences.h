/*
Original from Opengv.
Adaptation for the (central case) relative pose problem
*/


#ifndef GENRATECORRESPONDENCES_HPP_
#define GENRATECORRESPONDENCES_HPP_

#include <stdlib.h>
#include <vector>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <memory>

#include <Eigen/Core>

#include <opengv/math/cayley.hpp>

#include "EssentialTypes.h"
#include "CorrespondencesVectors.h"

namespace Essential
{

// Generate 3D point inside a frustum
Vector3 generateRandomPointTruncated(double FoV, double min_depth, double max_depth );
Vector3 generateRandomPointTruncated(double FoV, double min_depth, double max_depth, bool allow_coplanar, double max_X );

// add noise as explained in the paper(s)
Vector3 addNoise( double noiseLevel, Vector3 cleanPoint, double focal_length );
// generate random translation with the specified norm
Vector3 generateRandomTranslation( double max_parallax );
// generate random rotation with maxAngle
Matrix3 generateRandomRotation( double maxAngle );

Vector3 generateRandomPointNaive(double min_depth, double max_depth);


void createSyntheticExperiment(
    const size_t n_points,
    double noise,
    double outlier_fraction,
    double FoV,
    double max_parallax,
    double min_depth,
    double max_depth,
    Vector3 & translation,
    Matrix3 & rotation,
    bearing_vectors_t & points_correspondences,
    Eigen::MatrixXd & gt,
    std::vector<int> & indices_outliers,
    bool allow_coplanar,   // allow to generate synthetic scenes with coplanar points
    double max_X,  // this value is the absolute value. Max value for X-axis allowed 
    double min_epipolar_error_sq 
    );
    
    void createSyntheticExperiment(
    const size_t n_points,
    double noise,
    double outlier_fraction,
    double FoV,
    double max_parallax,
    double min_depth,
    double max_depth,
    Vector3 & translation,
    Matrix3 & rotation,
    bearing_vectors_t & points_correspondences,
    Eigen::MatrixXd & gt,
    std::vector<int> & indices_outliers,
    double min_epipolar_error_sq 
    );
    
    
    void createSyntheticExperiment(
    const size_t n_points,
    double noise,
    double outlier_fraction,
    double FoV,
    double max_parallax,
    double min_depth,
    double max_depth,
    Vector3 & translation,
    Matrix3 & rotation,
    bearing_vectors_t & points_correspondences,
    Eigen::MatrixXd & gt,
    std::vector<int> & indices_outliers,
    bool allow_coplanar,   // allow to generate synthetic scenes with coplanar points
    double max_X  // this value is the absolute value. Max value for X-axis allowed 
    );
    
    
void createSyntheticExperiment(
    const size_t n_points,
    double noise,
    double outlier_fraction,
    double FoV,
    double max_parallax,
    double min_depth,
    double max_depth,
    Vector3 & translation,
    Matrix3 & rotation,
    bearing_vectors_t & points_correspondences,
    Eigen::MatrixXd & gt,
    std::vector<int> & indices_outliers );


void createSyntheticExperiment(
    const size_t n_points,
    double noise,
    double outlier_fraction,
    double FoV,
    double parallax,
    double min_depth,
    double max_depth,
    Vector3 & translation,
    Matrix3 & rotation,
    bearing_vectors_t & points_correspondences,
    Eigen::MatrixXd & gt,
    std::vector<int> & indices_outliers,
    bool allow_coplanar,   // allow to generate synthetic scenes with coplanar points
    double max_X,  // this value is the absolute value. Max value for X-axis allowed
    double min_epipolar_error_sq,   // min epipolar error for the outliers
    double focal_length    // focal length in pixelss
    );
    
}  // end of essential namespace


#endif // end of GENRATECORRESPONDENCES_HPP_
