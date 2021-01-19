#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>  // for the file

// #include <opengv/types.hpp>
// #include <opengv/relative_pose/methods.hpp>
// #include <opengv/relative_pose/CentralRelativeAdapter.hpp>
#include "generateCorrespondences.h"

#include "Essential.h"
#include "EssentialTypes.h"
#include "EssentialUtils.h"




using namespace std;
using namespace Eigen;
using namespace opengv;
using namespace Essential;

int main(int argc, char** argv)
{
    std::cout << "Essential Matrix Estimation!\n";

  std::ofstream fout("init_rand.txt");
  double total_time = 0;
  double N_iter = 1;
  //set experiment parameters
  double noise = 5.5;
  size_t n_points = 200;
  double FoV = 100;  // in degrees
  double max_parallax = 2.;  // in meters
  double min_depth = 1.;     // in meters
  double max_depth = 8.;       // in meters
  double outlier_fraction = 0;


    std::srand(std::time(nullptr));

    for (int i = 0; i < N_iter; i++)
    {


       Vector3 translation;
       Matrix3 rotation;
       bearing_vectors_t points_correspondences;
       Eigen::MatrixXd points_3D(3, n_points);
       std::vector<int> indices_outliers(1, n_points);
       // FoV in degrees
       createSyntheticExperiment(n_points, noise, outlier_fraction, FoV,
                        max_parallax, min_depth, max_depth, translation, rotation,
                        points_correspondences, points_3D, indices_outliers);

      // n_point
      std::cout << "Computing essential matrix GT\n";
      Matrix3 E = computeEfromRt(rotation, translation);



      

      std::cout << "Ground truth rotation:\n" << rotation << std::endl;
      std::cout << "Ground truth translation Tgt:\n" << translation << std::endl;
      std::cout << " \n-------------------------\n";

      EssentialEstimationOptions options;
      options.chosen_initialisation = InitialisationMethod::PTS8_INIT;
      options.use_preconditioning = Preconditioner::Dominant_eigenvalues;
      options.verbose = 0;
      options.estimation_verbose=0;
      options.use_idx_relaxation = 0;
              
      EssentialClass my_essential_estimation(points_correspondences, options);

      // run the estimation
      
      std::cout << "Running algotihm on-manifold estimation\n";
      EssentialEstimationResult my_result = my_essential_estimation.getResults();
         

      my_essential_estimation.printResult(my_result);


      // Save results into a file
      /*
      {'iter','idx_iter','f_hat','d_hat','gradnorm','dual_gap','mu_min',
      'elapsed_init_time','elapsed_C_time','elapsed_8pt_time','elapsed_time_methods',
      'elapsed_iterative_time','elapsed_lagrange_time','elapsed_certifier_time','elapsed_estimation_time',
      'distT','distR','is_opt'}
      */
      fout << "iter: " << i;
      fout << " " << my_result.f_hat;
      fout << " " << my_result.d_hat;
      fout << " " << my_result.gradnorm;
      fout << " " << my_result.dual_gap;
      fout << " " << my_result.mu_min;
      fout << " " << my_result.elapsed_init_time;
      fout << " " << my_result.elapsed_C_time;
      fout << " " << my_result.elapsed_8pt_time;
      fout << " " << my_result.elapsed_time_methods;
      fout << " " << my_result.elapsed_iterative_time;
      fout << " " << my_result.elapsed_lagrange_time;
      fout << " " << my_result.elapsed_certifier_time;
      fout << " " << my_result.elapsed_estimation_time;
      fout << " " << distT(translation, my_result.t_opt);
      fout << " " << distR(rotation, my_result.R_opt);
      fout << " " << my_result.idx_relaxation; 
      
      int is_opt = -1;
      switch(my_result.certifier_status)
      {
        case EssentialEstimationStatus::RS_ITER_LIMIT:
            is_opt = 2;
            break;
        case EssentialEstimationStatus::GLOBAL_OPT:
            is_opt = 0;
            break;
        case EssentialEstimationStatus::NO_CERTIFIED:
            is_opt = 1;
            break;
        default:
            is_opt = -1;
      }


      fout << " " << is_opt << std::endl;

      total_time += my_result.elapsed_estimation_time;
      std::cout << "Total time: " << my_result.elapsed_estimation_time << std::endl;
      std::cout << "Error rotation: " << distR(rotation, my_result.R_opt) << std::endl;
      std::cout << "Error translation: " << distT(translation, my_result.t_opt) << std::endl;

  }
  // close file
  fout.close();
  std::cout << "Total time (total) [microsecs]: " << total_time << std::endl;
  std::cout << "Mean time (total) [microsecs]: " << total_time / N_iter << std::endl;

  return 0;

}
