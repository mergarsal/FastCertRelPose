#ifndef ESSENTIALMANIFOLD_H
#define ESSENTIALMANIFOLD_H

#include <random> // For sampling random points on the manifold

#include <Eigen/Dense>

#include "EssentialTypes.h"

/*Define the namespace*/
namespace Essential{

  class EssentialManifold{
  

      public:
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW 
        /* Default constructor */
        EssentialManifold(void){};

        /*Delete each component manifold*/
        ~EssentialManifold(void){};
        
        /// GEOMETRY
        /** Given a generic matrix A in R^{3 x 4}, this function computes the
        * projection of A onto R (closest point in the Frobenius norm sense).  */
        Matrix34 project(const Matrix34 &A) const;
        
        
        /** Given an element Y in M and a tangent vector V in T_Y(M), this function
       * computes the retraction along V at Y using the QR-based retraction
       * specified in eq. (4.8) of Absil et al.'s  "Optimization Algorithms on
       * Matrix Manifolds").
       */
      Matrix34 retract(const Matrix34 &Y, const Matrix34 &V) const;
            
      /* Projections for each manifold */
      /// Sphere
      Vector3 ProjSphere(const Vector3 &t, const Vector3 &Vt) const;
      /// Rotation  
      Matrix3 ProjRotation(const Matrix3& R, const Matrix3 & VR) const;
            
      // Product of the form: A * symm(B * C), where all the matrices are 3x3
      Matrix3 SymProduct(const Matrix3 & A, const Matrix3 B, const Matrix3 & C) const;
         
     /** Sample a random point on M, using the (optional) passed seed to initialize
       * the random number generator.  */
      Matrix34 random_sample(const std::default_random_engine::result_type &seed =
                               std::default_random_engine::default_seed) const;
                               
                               
                               /** Given an element Y in M and a matrix V in T_X(R^{p x kn}) (that is, a (p
   * x kn)-dimensional matrix V considered as an element of the tangent space to
   * the *entire* ambient Euclidean space at X), this function computes and
   * returns the projection of V onto T_X(M), the tangent space of M at X (cf.
   * eq. (42) in the SE-Sync tech report).*/
  Matrix34 Proj(const Matrix34 &Y, const Matrix34 &V) const;
    
      };
} /*end of Essential namespace*/
#endif // end of ESSENTIALMANIFOLD_H
