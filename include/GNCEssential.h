#pragma once

#include <chrono>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <experimental/optional>
#include <math.h>
#include <string>
#include <iomanip>
#include <algorithm>

// Essential estimation related
#include "EssentialTypes.h"
#include "EssentialUtils.h"
#include "EssentialProblem.h"
#include "EssentialManifold.h"
#include "EssentialVerification.h"
#include "Essential.h"

// generate correspondences
#include "CorrespondencesVectors.h"

// smooth optimization
#include "Optimization/Riemannian/Concepts.h"
#include "Optimization/Riemannian/TNT.h"
#include "Optimization/Riemannian/GradientDescent.h"
#include "Optimization/Convex/Concepts.h"
#include "Optimization/Convex/ProximalGradient.h"

// GNC related
#include "gncso/Smooth/gnc_smooth.h"

// gnc
using namespace gncso;
// Smooth
using namespace Optimization;
using namespace Optimization::Riemannian;

// essential matrix
namespace Essential {

  enum class GNCRobustFunction {
                  None = 0,
                  TLS,
                  GM,
                  TEMP,
                  WELSCH,
                  TUKEY,
                  Any,
              };


  struct GNCEssentialEstimationOptions : public GNCParams<double>, public EssentialEstimationOptions
  {
  GNCRobustFunction gnc_robust {GNCRobustFunction::Any};

  GNCEssentialEstimationOptions(): GNCParams<double>(), EssentialEstimationOptions() {}; 
  };


  struct GNCEssentialEstimationResult : public GNCResult<Matrix34, weights_t, double>, public EssentialEstimationResult
  {    
    // flag for valid estimation
    bool valid_estimation = false; 
    
    GNCEssentialEstimationResult(): EssentialEstimationResult(), GNCResult<Matrix34, weights_t, double>() {};
  };
  
  

  class GNCEssentialClass : public EssentialClass
  {
      public:
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW 
      /* Default constructor */
      GNCEssentialClass(void) : EssentialClass() {};

      GNCEssentialClass( const bearing_vectors_t & points,
                      const GNCEssentialEstimationOptions & options = GNCEssentialEstimationOptions(),
                      const Matrix3 & E_initial = Matrix3() ) :points_(points), E_initial_(E_initial), options_(options),
                      EssentialClass( points, options, E_initial) {};
                      
      GNCEssentialClass( Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors0,
                         Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors1,
                         Eigen::Matrix<double, 1, Eigen::Dynamic> & weights_observations,
                         const GNCEssentialEstimationOptions & options = GNCEssentialEstimationOptions(),
                         const Matrix3 & E_initial = Matrix3() ) : E_initial_(E_initial), options_(options),
                         EssentialClass( bearingVectors0, bearingVectors1, weights_observations, options, E_initial) 
                         {
                             bearing_vectors_t points = 
                             convertPointsVector32PointsOpt(bearingVectors0, bearingVectors1, weights_observations);
                             points_ = points;
                         }; 

      ~GNCEssentialClass(void) {};


      GNCEssentialEstimationResult getResultGNC(void);

      void printResultGNC(GNCEssentialEstimationResult & results);

      private:

         bearing_vectors_t points_;
         Matrix3 E_initial_;

         GNCEssentialEstimationOptions options_;


  };  //end of essential class}



}  // end of essential namespace
