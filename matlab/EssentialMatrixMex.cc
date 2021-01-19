/* Author:
    Mercedes Garcia Salguero

Code inspired by TEASER++ (see reference)

DISCL.:
THIS FILE IS NOT GNC

*/



#include <iostream>


#include <map>
#include "mex.h"
#include <Eigen/Core>

#include "EssentialMatrixMexUtils.h"

#include "Essential.h"
// #include "EssentialTypes.h"
// #include "EssentialUtils.h"


enum class INPUT_PARAMS : int
{
  obsi = 0,
  obsip,
  weights,
  init_method,
  precon,
  use_mult_certifiers,
  verbose
};

enum class OUTPUT_PARAMS : int
{
  R_est = 0,
  t_est,
  time_total,
  time_iterative,
  is_opt
};



typedef bool (*mexTypeCheckFunction)(const mxArray*);
const std::map<INPUT_PARAMS, mexTypeCheckFunction> INPUT_PARAMS_MAP{
    {INPUT_PARAMS::obsi, &is3NMatrix},
    {INPUT_PARAMS::obsip, &is3NMatrix},
    {INPUT_PARAMS::weights, &is1NMatrix},
    {INPUT_PARAMS::init_method, &isRealDoubleScalar},
    {INPUT_PARAMS::precon, &isRealDoubleScalar},
    {INPUT_PARAMS::use_mult_certifiers, &isRealDoubleScalar},
    {INPUT_PARAMS::verbose, &mxIsLogicalScalar},
};



const std::map<OUTPUT_PARAMS, mexTypeCheckFunction> OUTPUT_PARAMS_MAP{
    {OUTPUT_PARAMS::R_est, &isRealDoubleMatrix<3, 3>},
    {OUTPUT_PARAMS::t_est, &isRealDoubleMatrix<3, 1>},
    {OUTPUT_PARAMS::time_total, &isRealDoubleScalar},
    {OUTPUT_PARAMS::time_iterative, &isRealDoubleScalar},
    {OUTPUT_PARAMS::is_opt, &isRealDoubleScalar},
};




 void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {

   // TODO
  // Check for proper number of arguments
  
  if (nrhs != INPUT_PARAMS_MAP.size()) {
    mexErrMsgIdAndTxt("EssentialSolve:nargin", "Wrong number of input arguments.");
  }
  if (nlhs != OUTPUT_PARAMS_MAP.size()) {
    mexErrMsgIdAndTxt("EssentialSolve:nargin", "Wrong number of output arguments.");
  }
  
  

  
  // Check for proper input types
  for (const auto& pair : INPUT_PARAMS_MAP) {
    if (!pair.second(prhs[toUType(pair.first)])) {
      std::stringstream error_msg;
      error_msg << "Argument " << toUType(pair.first) + 1 << " has the wrong type.\n";
      mexErrMsgIdAndTxt("EssentialSolve:nargin", error_msg.str().c_str());
    }
  }
  

  mexPrintf("Arguments type checks passed.\n");
  mexEvalString("drawnow;");

  // Prepare parameters
  // Prepare source and destination Eigen point matrices
  Eigen::Matrix<double, 3, Eigen::Dynamic> obsi_eigen, obsip_eigen;
  Eigen::Matrix<double, 1, Eigen::Dynamic> weights_eigen;
  // conversion
  mex3NMatrixToEigenMatrix(prhs[toUType(INPUT_PARAMS::obsi)], &obsi_eigen);
  mex3NMatrixToEigenMatrix(prhs[toUType(INPUT_PARAMS::obsip)], &obsip_eigen);
  mex1NMatrixToEigenMatrix(prhs[toUType(INPUT_PARAMS::weights)], &weights_eigen);

  // check dimension
  assert(obsi_eigen.size() != 0);
  assert(obsip_eigen.size() != 0);
  assert(obsip_eigen.size() == obsi_eigen.size());
  assert(weights_eigen.size() != 0);

  // Other parameters
  // Note: init and precon have their own type
  auto init_method = static_cast<size_t>(*mxGetPr(prhs[toUType(INPUT_PARAMS::init_method)]));
  auto precon = static_cast<size_t>(*mxGetPr(prhs[toUType(INPUT_PARAMS::precon)]));
  auto use_mult_certifiers = static_cast<size_t>(*mxGetPr(prhs[toUType(INPUT_PARAMS::use_mult_certifiers)]));
  auto verbose = static_cast<bool>(*mxGetPr(prhs[toUType(INPUT_PARAMS::verbose)]));

  Essential::EssentialEstimationOptions options = Essential::EssentialEstimationOptions();
  // general verbose
  options.use_idx_relaxation = use_mult_certifiers;
  options.estimation_verbose = verbose;
  options.verbose = verbose;

  /* Initialization */
  /*  enum class InitialisationMethod {
                    PTS8_INIT = 0,
                    PTS7_INIT,                    
                    PTS5_INIT,
                    RANDOM_INIT,
                    EYE_INIT,
                    USER_SUPPLIED,
                };
  */
  options.chosen_initialisation = static_cast<Essential::InitialisationMethod>(init_method);

  /* Preconditioner */
     /*
          enum class Preconditioner {
                  None= 0,
                  Max_eigenvalue,
                  Dominant_eigenvalues, 
                  N,
                  Any,
               };
      */
  options.use_preconditioning = static_cast<Essential::Preconditioner>(precon);


  /**   SOLVER   **/
  // initialize estimator
  Essential::EssentialClass my_essential_estimation(obsi_eigen, obsip_eigen, weights_eigen, options);


  mexPrintf("Start essential matrix estimator.\n");
  mexEvalString("drawnow;");

  // Solve
  Essential::EssentialEstimationResult my_result = my_essential_estimation.getResults();
  /*
   end of SOLVER
  */

  // Populate output E matrix
  plhs[toUType(OUTPUT_PARAMS::R_est)] = mxCreateDoubleMatrix(3, 3, mxREAL);
  Eigen::Map<Eigen::Matrix3d> R_map(mxGetPr(plhs[toUType(OUTPUT_PARAMS::R_est)]), 3, 3);
  R_map = my_result.R_opt;

  // Populate output T vector
  plhs[toUType(OUTPUT_PARAMS::t_est)] = mxCreateDoubleMatrix(3, 1, mxREAL);
  Eigen::Map<Eigen::Matrix<double, 3, 1>> t_map(mxGetPr(plhs[toUType(OUTPUT_PARAMS::t_est)]), 3, 1);
  t_map = my_result.t_opt;


  // Populate outputs
  // iterative algorithm
  plhs[toUType(OUTPUT_PARAMS::time_iterative)] = mxCreateDoubleScalar(my_result.elapsed_iterative_time);
  plhs[toUType(OUTPUT_PARAMS::time_total)] = mxCreateDoubleScalar(my_result.elapsed_estimation_time);

  // is opt?
  double is_opt_double = -1;
  switch(my_result.certifier_status)
  {
    case (Essential::EssentialEstimationStatus::GLOBAL_OPT):
    {
      mexPrintf("Solution is optimal!.\n");
      is_opt_double = 1;
      break;
    }
    case (Essential::EssentialEstimationStatus::NO_CERTIFIED):
    {
      mexPrintf("We could not certify the solution.\n");
      is_opt_double = 0;
      break;
    }
    default:
    {
      mexPrintf("Max. nr iterations in solver. Try increasing the parameter.\n");
      break;
    }
  }

  plhs[toUType(OUTPUT_PARAMS::is_opt)] = mxCreateDoubleScalar(is_opt_double);
}
