#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

#include "generateCorrespondences.h"
#include "Essential.h"
#include "EssentialTypes.h"
#include "EssentialUtils.h"

#include "GNCEssential.h"


using namespace std;
using namespace Eigen;
using namespace Essential;



int main(int argc, char** argv)
{
    std::cout << "Essential Matrix Estimation with GNC GM!\n";

  // save general results
  std::ofstream fout("GNC_TNT_res.txt");
  // save weights from the GNC
  std::ofstream fout_w("GNC_TNT_weights.txt");
  
  double N_iter = 1;  
  //set experiment parameters
  double noise = 0.1;
  size_t n_points = 100;
  double FoV = 150;  // in degrees
  double parallax = 2.0;  // in meters
  double min_depth = 1.;     // in meters
  double max_depth = 8.;       // in meters
  double outlier_fraction = 0.10;


   std::srand(std::time(nullptr));

    for (int i = 0; i < N_iter; i++)
    {
          Vector3 translation;
          Matrix3 rotation;
          bearing_vectors_t points_correspondences;
          Eigen::MatrixXd points_3D(3, n_points);
          std::vector<int> indices_outliers(1, n_points);
           createSyntheticExperiment(n_points, noise, outlier_fraction, FoV, parallax, min_depth, max_depth, translation, rotation,
                            points_correspondences, points_3D, indices_outliers, false, 20, 0.3 );
                            
    
        
          // n_point
          Matrix3 E = computeEfromRt(rotation, translation);


      
      GNCEssentialEstimationOptions options = GNCEssentialEstimationOptions();
      /* For GNC */ 
      options.gnc_robust = GNCRobustFunction::WELSCH;
      options.GNC_verbose = 0;
      options.gnc_factor = 1.10;       // for GM
      options.cost_diff_threshold = 0.000001;  // stop criteria. if the diff between two cost is lower than this value, stop
      options.max_res_tol_sq = 0.00001;  
      options.max_inner_iterations = 2;  // with 4 (even less) you have enough
      
      // options.mu_threshold = 1+1e-08;           // for GM
      // options.mu_threshold = 0.001;             // for TLS
      options.GNClog_iterates = true;              // log weights. NOTE: deactive this because memory
      options.inliers_threshold = 0.9;             // for GM
      options.nr_min_points = 12;
            
      
      /* For the TNT solver */
      
      options.chosen_initialisation = InitialisationMethod::PTS8_INIT;       
      options.use_preconditioning = Preconditioner::Any;
      options.estimation_verbose = 0;    
      
      
      
      // Instance of GNC estimator
      GNCEssentialClass my_essential_estimation(points_correspondences, options, E);

      // run the estimation
      GNCEssentialEstimationResult my_result = my_essential_estimation.getResultGNC();
      std::cout << "GT R:\n" << rotation << std::endl;
      std::cout << "Estimated R:\n" << my_result.R_opt << std::endl;
      std::cout << "Solution is valid?: " << my_result.valid_estimation << std::endl;

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


      fout << " " << is_opt;
      
      fout << " " << my_result.valid_estimation << std::endl;
      

      std::cout << "Total time: " << my_result.elapsed_estimation_time << std::endl;    
      std::cout << "Error rotation: " << distR(rotation, my_result.R_opt) << std::endl;
      std::cout << "Error translation: " << distT(translation, my_result.t_opt) << std::endl;
      std::cout << "GNC time: " << my_result.elapsed_iterative_time << std::endl;
      
      
      
      /* Saving weights */ 
      // handle: fout_w
      // format: one file per iteration
      // number of outer iterations: my_result.nr_outer_iterations
      for (size_t id_outer = 0; id_outer < my_result.nr_outer_iterations; id_outer++)
      {
        // extract weights    
        fout_w << id_outer;
        for (size_t id_w = 0; id_w < my_result.intermediate_outer_results[id_outer].weights.size(); id_w ++)
        {
                fout_w << "," << my_result.intermediate_outer_results[id_outer].weights[id_w];
        }
        // add f(hat(x))
        fout_w << "," << my_result.intermediate_outer_results[id_outer].f;
        fout_w << std::endl;
      }
      
     
      

  }
  // close file for general results
  fout.close(); 
  // close file for weights
  fout_w.close();  

  return 0;

}
