#ifndef _ESSENTIAL_H_
#define _ESSENTIAL_H_

#include <functional>

#include <chrono>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <math.h>

#include "EssentialTypes.h"
#include "EssentialUtils.h"
#include "EssentialProblem.h"
#include "EssentialManifold.h"
#include "EssentialVerification.h"

#include "CorrespondencesVectors.h"
#include <experimental/optional>
#include <Optimization/Riemannian/TNT.h>
#include "Optimization/Convex/Concepts.h"

#include "Optimization/Riemannian/GradientDescent.h"
#include "Optimization/Convex/ProximalGradient.h"


#include <string>
#include <iomanip>
#include <algorithm>


namespace Essential
{
    enum class InitialisationMethod {
                    PTS8_INIT = 0,
                    PTS7_INIT,
                    PTS5_INIT,
                    RANDOM_INIT,
                    EYE_INIT,
                    USER_SUPPLIED,
                };


    typedef Optimization::Riemannian::TNTUserFunction<Matrix34, Matrix34, double, Matrix34>
                        EssentialTNTUserFunction;

        /** This struct contains the various parameters that control the algorithm 
        Based on: SESync struct*/
    struct EssentialEstimationOptions {

      /// OPTIMIZATION STOPPING CRITERIA
      /** Stopping tolerance for the norm of the Riemannian gradient */
      double tol_grad_norm = 1.e-9;

       /** Stopping criterion based upon the norm of an accepted update step */
      double preconditioned_grad_norm_tol = 1.e-9;

      /** Stopping criterion based upon the relative decrease in function value */
      double rel_func_decrease_tol = 1e-10;

       /** Stopping criterion based upon the norm of an accepted update step */
      double stepsize_tol = 1e-03;

       /** Gradient tolerance for the truncated preconditioned conjugate gradient
   * solver: stop if ||g|| < kappa * ||g_0||.  This parameter should be in the
   * range (0,1). */
      double STPCG_kappa = 0.7;

      /** Gradient tolerance based upon a fractional-power reduction in the norm of
   * the gradient: stop if ||g|| < ||kappa||^{1+ theta}.  This value should be
   * positive, and controls the asymptotic convergence rate of the
   * truncated-Newton trust-region solver: specifically, for theta > 0, the TNT
   * algorithm converges q-superlinearly with order (1+theta). */
      double STPCG_theta = .5;

      /** Maximum permitted number of (outer) iterations of the RTR algorithm */
      unsigned int max_RTR_iterations = 2000;
      /** Maximum number of inner (truncated conjugate-gradient) iterations to
      * perform per out iteration */
      unsigned int max_tCG_iterations = 5;

      /// ESSENTIAL ESTIMATION PARAMS.
      // Absolute tolerance for the dual gap f - d
      double dual_gap_tol = 1e-10;

      // tol for the minimum eigenvalue of the (penalised) matrix
      double eig_min_tol = -0.02;

      // variable for the type of initialisation: 0 -> 8pts, 1->random, 2->eye, 3-> user supplied
      InitialisationMethod chosen_initialisation = InitialisationMethod::PTS8_INIT;

      /** Whether to use preconditioning in the trust regions solver */
      Preconditioner use_preconditioning = Preconditioner::Any;

      /* Whether to try *all* the relaxations till one comes back positive 
      or we run out of relaxations without a positive certifier */
      size_t use_idx_relaxation = 0;  // Note: use_idx_relaxation \in {0, 5}
      
      bool use_all_relaxations = true;  // when FALSE, it  runs relaxation 'use_idx_relaxation'
                                 // when TRUE, it runs UP TO relaxation 'use_idx_relaxation'

      // verbose = {0, 1}
      unsigned int estimation_verbose = 0;
      
      unsigned int verbose = 0;  // for TNT


      /// DEFAULT CONSTRUCTOR with default values

       EssentialEstimationOptions() {};
       

};  // end of struct: EssentialEstimationOptions




enum class EssentialEstimationStatus {
  /** The algorithm converged to a certified global optimum */
  GLOBAL_OPT = 0,

   // The primal point couldnt be certified
  NO_CERTIFIED,

  /** The algorithm exhausted the maximum number of iterations of the Riemannian
    optimization*/
  RS_ITER_LIMIT
};




/** This struct contains the output of the Essential Estimation */
struct EssentialEstimationResult {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW 
  
  // Primal objective value
  double f_hat;

  // Dual objective value
  double d_hat;

  // dual gap
  double dual_gap;
  // The norm of the Riemannian gradient at x_hat
  double gradnorm;

  // The minimum eigenvalue of the penalised matrix M
  double mu_min;

  /** The corresponding lagrange multipliers */
  Vector6 lagrange_multipliers;

  /* condition number k(hessian) */
  double condition_number;

   // Elapsed time for the initialisation
   double elapsed_init_time;
   // Elapsed time for C
   double elapsed_C_time;
   // Elapsed time for computing the 8 point init + preconditioner
   double elapsed_8pt_time;
   // Elapsed time for computing the initial estimation
   double elapsed_time_methods;
   // Elapsed time for the optimization on manifold
   double elapsed_iterative_time;
   // Elapsed time: lagrange multipliers
   double elapsed_lagrange_time;
   // Elapsed time for the certification
   double elapsed_certifier_time;
   // Elapsed time for the whole estimation: initialisation + optimization + verification
   double elapsed_estimation_time;

   // index of the relaxation which could certify optimality. 
   // idx_relaxation \in {0, 5}; idx_relaxation = 10 (no certified)
   size_t idx_relaxation; 
   
   // Output rotation matrix
   Matrix3 E_opt;

   // Output rotation matrix
   Matrix3 R_opt;
   // Output translation vector
   Vector3 t_opt;
   // Status
   EssentialEstimationStatus certifier_status;

    /* Default constructor */
   EssentialEstimationResult() {};

}; // end of EssentialEstimationResult struct




class EssentialClass
{
    public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW 
    /* Default constructor */
    EssentialClass(void){};

    EssentialClass( const bearing_vectors_t & points,
                    const EssentialEstimationOptions & options = EssentialEstimationOptions(),
                    const Matrix3 & E_initial = Matrix3()):points_(points), E_initial_(E_initial), options_(options) {};
                    
    EssentialClass( Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors0,
                    Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors1,
                    Eigen::Matrix<double, 1, Eigen::Dynamic> & weights_observations,
                    const EssentialEstimationOptions & options = EssentialEstimationOptions(),
                    const Matrix3 & E_initial = Matrix3()): E_initial_(E_initial), options_(options) 
                    {
                        bearing_vectors_t points = 
                                       convertPointsVector32PointsOpt(bearingVectors0, bearingVectors1, weights_observations);
                        points_ = points;
                    };        
                           

    ~EssentialClass(void){};

    EssentialEstimationResult getResults(void);

    void printResult(EssentialEstimationResult & results);

    private:
        bearing_vectors_t points_;
        Matrix3 E_initial_;
        EssentialEstimationOptions options_;

};  //end of essential class

}  // end of essential namespace
#endif // _ESSENTIAL_H_
