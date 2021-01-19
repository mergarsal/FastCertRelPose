#pragma once


#include <Eigen/SVD> // required for SVD decomposition in unlift method
#include <Eigen/Core>
#include <vector>

#include <string>

#include "EssentialTypes.h"
#include "CorrespondencesVectors.h"

#include <opengv/types.hpp>  // Not sure if we need this
#include <opengv/relative_pose/CentralRelativeAdapter.hpp>
#include <opengv/relative_pose/methods.hpp>

#define SQRT2 1.41421356237
#define R180PI 57.2957795131

/* Define namespace */
namespace Essential{



                // Return the cross-product for 3D vectors as a 3x3 matrix
                Matrix3 cross(const Vector3 & t);

                // Return the 9D vector for a 3x3 matrix (column-wise)
                Vector9 vec(Matrix3 & M);

                // Return the 9D vector for a 3x3 matrix (column-wise)
                Eigen::MatrixXf skew(const Eigen::MatrixXf & M);

        
                // only points. Basic
                Matrix3 initialize8pts(const bearing_vectors_t & points);

                // 8pts initialization
                Matrix3 initialize8pts(const bearing_vectors_t & points, const Matrix9 & C, Matrix3 & Jacobi_precon, Preconditioner use_precon);

                // Same as above but with C instead of the points.
                // In practice, they are the same
                Matrix3 initialize8pts(const Matrix9 & C , Matrix3 & Jacobi_precon, Preconditioner use_precon);

                /* minimal solvers */ 
                Matrix3 minimalinitialize8pts(const bearing_vectors_t & points, bool use_all_corr = false);
                // Use all the points 
                Matrix3 minimalinitialize8ptsAll(const bearing_vectors_t & points);
                // use only 8 points 
                Matrix3 minimalinitialize8ptsMinimal(const bearing_vectors_t & points);
                
                
                // 5pts initialization
                Matrix3 initialize5pts(const bearing_vectors_t & points, bool use_all_corr = false);
                
                // use all points 
                Matrix3 initialize5ptsAll(const bearing_vectors_t & points); 
                
                // use only 5 
                Matrix3 initialize5ptsMinimal(const bearing_vectors_t & points);
                

                // triangulate points 
                Vector3 triangulatePoint(const Matrix3 & R, const Vector3 & t, const Vector3 & f0, const Vector3& f1);
                
                
                // 7pts initialization
                Matrix3 initialize7pts(const bearing_vectors_t & points, bool use_all_corr = false); 
                
                // use all points 
                Matrix3 initialize7ptsAll(const bearing_vectors_t & points); 
                
                // use only 7 points
                Matrix3 initialize7ptsMinimal(const bearing_vectors_t & points);
                

                // Project to rotation group 
                Matrix3 projectToRotation(const Matrix3 & R_hat);
                
                // Project to Me
                Matrix3 projectToEssentialManifold(const Matrix3 & E_hat);

                // Construct That such that: vec(E) = e = Rhat * t
                Matrix39 createRhat(const Matrix3& R);

                // Construct That such that: vec(E) = e = that * r = that * vec(R)
                Matrix9 createThat(const Vector3 & t);

                // construct Mt, such that: f(E) = r^T * Mrt * r, r = vec(R), Mt = t_hat^T * C * t_hat
                Matrix3 createMatrixT(const Matrix9 & C, const Matrix3 & R);

                // construct Mr, such that: f(E) = t^T * Mr * t, Mr = R_hat^T * C * R_hat
                Matrix9 createMatrixR(const Matrix9 & C, const Vector3 & t);

                // construct C = sum_{i=1}^N C_i
                Matrix9 constructDataMatrix(const bearing_vectors_t & bearing_vectors);

                // construct C = sum_{i=1}^N C_i
                Matrix9 constructDataMatrix(const bearing_vectors_t & bearing_vectors, const weights_t & new_weights);


                // Construct B's matrices, use in the derivate
                void constructBMatricesForMixedDiftR(Matrix9& B1, Matrix9 & B2, Matrix9 & B3);

                // Construct matrix for mixed derivate
                Matrix39 constructMixedDifTtr(const Matrix9 & B1, const Matrix9 & B2,
                                const Matrix9 & B3,const Matrix9 & C , const Vector3 & t, Matrix3 & R);


                // Symmetrize the given matrix A = matrix_data. output: 0.5*(A + A^T)
                Eigen::MatrixXf symmetrize(const Eigen::MatrixXf & matrix_data);


                // Extract R, t from E:
                void computeRtfromE(const bearing_vectors_t & points, const Matrix3& E_hat, Matrix3 & R, Vector3 & t );

                // Compute E from R, t
                Matrix3 computeEfromRt(const Matrix3 & R, const Vector3 & t);

                Matrix3 computeEfromRt(const Matrix34& Rt);

                // Transform from Vector3 (eigen) to our struct
                bearing_vectors_t convertPointsVector32PointsOpt(std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> & bearingVectors0,
                                                                 std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> & bearingVectors1);

                bearing_vectors_t convertPointsVector32PointsOpt(std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> & bearingVectors0,
                                                                std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> & bearingVectors1, 
                                                                std::vector<double> & weights_observations);

                bearing_vectors_t convertPointsVector32PointsOpt(Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors0,
                                                                 Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors1,
                                                                 Eigen::Matrix<double, 1, Eigen::Dynamic> & weights_observations);

                // compute distance (in geodesic sense) between two rotations (See paper) [degrees]
                double distR(const Matrix3 & R_gt, const Matrix3 & R_est);
                // Compute distance (in direction) between two (unitary) translation vectors [degrees]
                double distT(const Vector3 & T_gt, const Vector3 & T_est);
                // Compute similarity between two essential matrices [degrees]
                double distE(const Matrix3 & E_gt, const Matrix3 & E_est);
                // Evaluate essential matrix, rotation and translation direction (everything in [degreees])

                std::vector<int> selectKPoints(int max_n, int k); 
                
                void evaluateRTE(const Matrix3 & R_gt, const Matrix3 & R_est,
                            const Vector3 & T_gt, const Vector3 & T_est,
                            const Matrix3 & E_gt, const Matrix3 & E_est,
                            double& dist_R, double & dist_T, double & dist_E);
} // end of essential namespace


