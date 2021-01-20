#pragma once



/* Essential manifold related */
#include "EssentialTypes.h"
#include "EssentialManifold.h"
#include "CorrespondencesVectors.h"
#include "EssentialUtils.h"

#include <Eigen/Dense>

namespace Essential {


struct ProblemCachedMatrices{
EIGEN_MAKE_ALIGNED_OPERATOR_NEW 
        Matrix34 NablaF_Y;
        Matrix3 Mt;
        Matrix9 Mr;
        /// DEFAULT CONSTRUCTOR with default values
        ProblemCachedMatrices(  const Matrix34 & Nabla = Matrix34::Identity(),
                                const Matrix3 &  mt = Matrix3::Identity(),
                                const Matrix9 & mr = Matrix9::Identity()) :
                                NablaF_Y(Nabla), Mt(mt), Mr(mr){}
};


class EssentialProblem{
public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW 
        EssentialProblem(){};  // Default

        // Constructor using two vectors of 3XN corresponding features
        EssentialProblem(const bearing_vectors_t & bearing_vectors, Preconditioner use_precon = Preconditioner::None);

        void setNumberPoints(const bearing_vectors_t& bearing_vectors);

        EssentialProblem(const Matrix9 & data_matrix_C, Preconditioner use_precon = Preconditioner::None);

        ~EssentialProblem();

         Matrix9 getDataMatrixC(void) {return data_matrix_C_;}
         
         void setPointCorrespondences(const bearing_vectors_t & bearing_vectors) {point_correspondences_ = bearing_vectors; number_points_ = bearing_vectors.size();}


         void setMatrixPrecon(Matrix3 & matrix_precon) {Matrix_precon_ = matrix_precon;}


         // Initialize solver with 8 points
         Matrix34 initializeSolver8pts(void);

         // This preconditioner uses P = Mt
         Matrix3 invMtPrecon(Matrix3 & R_init);

         // Pseudo jacobi preconditioner based on the three largest eigenvalues of C
         Matrix3 computePseudoJacobiPrecon(void);


        /// ACCESSORS

        /** Get a const pointer to the SO(3) x S(2) product manifold */
         const EssentialManifold& getEssentialManifold() const {
            return domain_;
          }


          Matrix9 getDataMatrix(void) {return data_matrix_C_;}

          void setMr(const Matrix9 & M) {Mr_ = M;};
          void setMt(const Matrix3 & M) {Mt_ = M;};

          Matrix3 getMt(void) {return Mt_; };
          Matrix9 getMr(void) {return Mr_; };
           /// OPTIMIZATION AND GEOMETRY

          /** Given a matrix Y, this function computes and returns F(Y), the value of
   * the objective evaluated at X */
  double evaluate_objective(const Matrix34 &Y) const;

  double evaluate_objective(const Matrix34 &Y, ProblemCachedMatrices & problem_matrices) const;

    /** Given a matrix Y, this function computes and returns nabla F(Y), the
   * *Euclidean* gradient of F at Y. */
  Matrix34 Euclidean_gradient(const Matrix34 &Y) const;

  Matrix34 Euclidean_gradient(const Matrix34 &Y, const ProblemCachedMatrices & problem_matrices) const;

  /** Given a matrix Y in the domain D of the SE-Sync optimization problem and
   * the *Euclidean* gradient nabla F(Y) at Y, this function computes and
   * returns the *Riemannian* gradient grad F(Y) of F at Y */

   Matrix34 Riemannian_gradient(const Matrix34 &Y, const ProblemCachedMatrices & problem_matrices) const;

   Matrix34 Riemannian_gradient(const Matrix34 &Y, const Matrix34 &nablaF_Y) const;

   Matrix34 Riemannian_gradient(const Matrix34 &Y) const;

  /* Preconditioner */
  Matrix34 precondition(const Matrix34& X, const Matrix34 & Xdot) const;

  Matrix34 preconditionRight(const Matrix34& X, const Matrix34 & Xdot) const;

  /** Given a matrix Y in the domain D of the SE-Sync optimization problem, the
   * *Euclidean* gradient nablaF_Y of F at Y, and a tangent vector dotY in
   * T_D(Y), the tangent space of the domain of the optimization problem at Y,
   * this function computes and returns Hess F(Y)[dotY], the action of the
   * Riemannian Hessian on dotY */

   Matrix34 Riemannian_Hessian_vector_product(const Matrix34 &Y,
                                                   const ProblemCachedMatrices & problem_matrices,
                                                   const Matrix34 &dotY) const;

   Matrix34 Riemannian_Hessian_vector_product(const Matrix34 &Y,
                                           const Matrix34 &nablaF_Y,
                                           const Matrix34 &dotY) const;


   Matrix34 Riemannian_Hessian_vector_product(const Matrix34 &Y,
                                            const Matrix34 &dotY) const;

   Matrix34 Riemannian_Hessian_vector_product(const Matrix34 &Y,
               const Matrix3 & Mt,
               const Matrix9 & Mr,
               const Matrix34 &dotY) const;

    /** Given a matrix Y in the domain D of the SE-Sync optimization problem and a
  tangent vector dotY in T_Y(E), the tangent space of Y considered as a generic
  matrix, this function computes and returns the orthogonal projection of dotY
  onto T_D(Y), the tangent space of the domain D at Y*/
  Matrix34 tangent_space_projection(const Matrix34 &Y, const Matrix34 &dotY) const;

  /** Given a matrix Y in the domain D of the SE-Sync optimization problem and a
   * tangent vector dotY in T_D(Y), this function returns the point Yplus in D
   * obtained by retracting along dotY */
  Matrix34 retract(const Matrix34 &Y, const Matrix34 &dotY) const;

  Matrix34 random_sample() const;

  // compute Residuals
  weights_t computeResiduals(const Matrix34 & Rt) const;

  // update points with weights
  void updateWeights(const weights_t & new_weights);  // This function modifies the weights

private:

  size_t number_points_;

  bearing_vectors_t point_correspondences_;  // with weights!!

  Preconditioner use_precon_;

  /** The product manifold of SO(3) X S(2) ~ E(3) that is the domain of our method */
  EssentialManifold domain_;

  Matrix9 data_matrix_C_;

  Matrix3 Mt_;

  Matrix9 Mr_;

  Matrix3 Matrix_precon_;

  // Matrices for the mixed term in the hessian
  Matrix9 B1_, B2_, B3_;

  // matrices for the preconditioner: B^T * C * B
  Matrix93 B_;

  Matrix4 Matrix_precon_right_;

}; // end of Essential problem class




}  // end of Essential namespace


