
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

// Eigen
#include <Eigen/Dense>

// Essential
#include "generateCorrespondences.h"
#include "Essential.h"
#include "EssentialTypes.h"
#include "EssentialUtils.h"

// GNC
#include "gncso/Base/GNC.h"

#define N_POINTS 100
typedef double Scalar;
typedef Eigen::Matrix<Scalar, 3, N_POINTS> Points3D;
typedef Eigen::Matrix<Scalar, 1, N_POINTS> VectorN;
typedef Eigen::Matrix<Scalar, 1, N_POINTS> Weights;



using namespace std;
using namespace Eigen;
using namespace Essential;

int main(int argc, char** argv)
{

  std::cout << "Naive Test Essential Matrix\n";

  double total_time = 0;
  double N_iter = 1;
  //set experiment parameters
  double noise = 0.;
  double FoV = 100;  // in degrees
  double max_parallax = 2.;  // in meters
  double min_depth = 1.;     // in meters
  double max_depth = 8.;       // in meters
  double outlier_fraction = 0.5;

  // random seed
  std::srand(std::time(nullptr));

    for (int i = 0; i < N_iter; i++)
    {
          Vector3 translation;
          Matrix3 rotation;
          bearing_vectors_t points_str;
          Eigen::MatrixXd points_3D(3, N_POINTS);
          std::vector<int> indices_outliers(1, N_POINTS);
          createSyntheticExperiment(N_POINTS, noise, outlier_fraction, FoV, max_parallax, min_depth, max_depth, translation, rotation, points_str, points_3D, indices_outliers, 0.02 );


          Matrix3 E_est, E_ref = computeEfromRt(rotation, translation);

       // Define the estimation function
       baseOpt::VariableEstimation<Matrix3, Weights, bearing_vectors_t> compute_svd_fcn = [](const Matrix3& X, const Weights& weights, const bearing_vectors_t& points)
      {
          std::cout << "weight: " << weights << std::endl;
          Matrix9 C = constructDataMatrix(points, weights);

          Eigen::JacobiSVD<Matrix9> svd_C(C, Eigen::ComputeFullU | Eigen::ComputeFullV);

          // force the solution to be an essential matrix
          Vector9 e_null = svd_C.matrixV().col(8);
          Matrix3 E_null = Eigen::Map<Matrix3>(e_null.data(), 3, 3);
          Eigen::JacobiSVD<Matrix3>svd(E_null, Eigen::ComputeFullU | Eigen::ComputeFullV);
          Eigen::DiagonalMatrix<double, 3> d(1, 1, 0);
          std::cout << "\nEiegenvalues: " << svd_C.singularValues() << std::endl << std::endl;

          Matrix3 U = svd.matrixU();
          Matrix3 V = svd.matrixV();

          // force the solution to have SVD with rotation matrices
          if (U.determinant() < 0)    U.col(2) *= -1;
          if (V.determinant() < 0)    V.col(2) *= -1;

          return (U * d * (V).transpose());
       };

       // Define the cost function

       baseOpt::ComputeCost<Matrix3, Weights, Scalar, bearing_vectors_t> compute_cost_fcn = [](const Matrix3& X, const Weights& weights, const Weights& residuals_sq, const bearing_vectors_t& points )
       {
        Scalar cost = 0;
        for (size_t i = 0; i < N_POINTS; i ++)    cost += weights(i) * residuals_sq(i);

        return (cost);
      };

       // Define how to compute the residuals ^2
       baseOpt::ComputeResiduals<Matrix3, Weights, bearing_vectors_t> compute_residuals_fcn = [](const Matrix3& X, const bearing_vectors_t& points)
       {
         Weights residual;

         for (size_t i = 0; i < N_POINTS; i++)
         {
           Vector3 v0 = points[i].bearing_vector_0;
           Vector3 v1 = points[i].bearing_vector_1;
           residual(i) = (pow((v0.transpose() * X * v1), 2));
         }

         return (residual);
       };



       // Call solver w/t noise
       gncso::GNCParams<Scalar> options = gncso::GNCParams<Scalar>();
         options.inliers_threshold = 0.8;       // threshold to apply to the inliers. IF weight(i) > inliers_threshold, then 'i' is an inlier
         options.cost_diff_threshold = 0.0001;  // stop criteria. if the diff between two cost is lower than this value, stop
         options.max_res_tol_sq = 0.02 * 0.02;  // maximum tolerance allowed (square)
         options.mu_threshold = 1+1e-08;          // for GM
         options.GNC_verbose = 1;
         // Initial estimations

         Weights weights_initial;
         weights_initial.setOnes(1, N_POINTS);
        
         

         std::cout << compute_residuals_fcn(E_ref, points_str) << std::endl;
         Matrix3 E_init = compute_svd_fcn(Matrix3::Identity(), weights_initial, points_str);
         
                
         
         

   std::cout << "Solving problem without outliers nor noise...\n";
   // Solve the problem!!!
   gncso::GNCResult<Matrix3, Weights, Scalar> results = gncso::GMGNC<Matrix3, Weights, Scalar, bearing_vectors_t>(compute_svd_fcn, compute_cost_fcn, compute_residuals_fcn, E_init, weights_initial, points_str, std::experimental::nullopt, std::experimental::nullopt, options);
   // extract solution
   E_est = results.x;
   std::cout << "Result without outliers\n------------\n";
   std::cout << "Ground truth rotation matrix:\n" << E_ref << std::endl;
   std::cout << "Estimated rotation matrix:\n" << E_est << std::endl;
   std::cout << "Geodesic distance between both rotations: " << distE(E_ref, E_est) << std::endl;
   std::cout << "Detected inliers: " << results.set_inliers << std::endl;

  }  // end for each iteration

  return 0;

}
