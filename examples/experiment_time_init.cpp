
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <opengv/relative_pose/methods.hpp>
#include <opengv/relative_pose/CentralRelativeAdapter.hpp>
#include <sstream>
#include <fstream>

#include <opengv/types.hpp>

#include "Essential.h"
#include "generateCorrespondences.h"
#include "EssentialTypes.h"
#include "EssentialUtils.h"


#include <fstream>  // for the file

using namespace std;
using namespace Eigen;
using namespace opengv;
using namespace Essential;

int main(int argc, char** argv)
{
  std::cout << "Essential Matrix Estimation!\n";


  double total_time = 0;
  double N_iter = 100;
  //set experiment parameters
  double noise = 0.5;
  double FoV = 100;  // in degrees
  double max_parallax = 2.;  // in meters
  double min_depth = 1.0;     // in meters
  double max_depth = 8.;       // in meters
  double outlier_fraction = 0;
  size_t n_noise = 4;
  size_t n_points_id = 12;
  Eigen::MatrixXd array_n_points(20, 1);

  array_n_points << 8, 9, 10, 11, 12, 13, 14, 15, 20, 40, 100, 200;


  std::srand(std::time(nullptr));
  size_t n_init = 3;

  for (size_t i = 0; i < n_init ; i++)
  {
    InitialisationMethod init_method = static_cast<InitialisationMethod>(i);

    for (size_t points_id=0; points_id < n_points_id; points_id++)
    {
        size_t n_points = array_n_points(points_id);

        auto name_file = "test_noise_05_N" + std::to_string(points_id) + "_init_" + std::to_string(i) + ".txt";

        std::ofstream fout(name_file);

        for (int i = 0; i < N_iter; i++)
        {


            Vector3 translation;
                  Matrix3 rotation;
                  bearing_vectors_t points_correspondences;
                  Eigen::MatrixXd points_3D(3, n_points);
                  std::vector<int> indices_outliers(1, n_points);

               createSyntheticExperiment(n_points, noise, outlier_fraction, FoV,
                                max_parallax, min_depth, max_depth, translation, rotation,
                                points_correspondences, points_3D, indices_outliers );

              EssentialEstimationOptions options;
              options.chosen_initialisation = init_method;
              EssentialClass my_essential_estimation(points_correspondences, options);


              EssentialEstimationResult my_result = my_essential_estimation.getResults();

       

              // Save results into a file

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

       }
      // close file
      fout.close();

    }  // end for each point level

}  // end for each init
  return 0;

}
