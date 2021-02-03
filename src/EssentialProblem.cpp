#include "EssentialProblem.h"


namespace Essential{
        EssentialProblem::EssentialProblem(const bearing_vectors_t& bearing_vectors, Preconditioner use_precon) : use_precon_(use_precon)
        {

                point_correspondences_ = bearing_vectors;
                data_matrix_C_ = constructDataMatrix(bearing_vectors);
                number_points_ = bearing_vectors.size();

                Mt_.setZero();
                Mr_.setZero();
                constructBMatricesForMixedDiftR(B1_, B2_, B3_);


        }; //end of constructor
        // Delete this

        EssentialProblem::EssentialProblem(const Matrix9& data_C, Preconditioner use_precon) : use_precon_(use_precon)
        {
                data_matrix_C_ = data_C;

                Mt_.setZero();
                Mr_.setZero();
                constructBMatricesForMixedDiftR(B1_, B2_, B3_);

        }; //end of constructor

        void EssentialProblem::setNumberPoints(const bearing_vectors_t& bearing_vectors)
        {
          number_points_ = bearing_vectors.size();
        }

        EssentialProblem::~EssentialProblem(void)
        { }
        

        weights_t EssentialProblem::computeResiduals(const Matrix34 & Rt) const
        {
            weights_t residual;
            residual.resize(number_points_);

            Matrix3 E = computeEfromRt(Rt);

            for(size_t i = 0; i < number_points_; i++)
            {
                    Vector3 v1 = point_correspondences_[i].bearing_vector_1;
                    Vector3 v0 = point_correspondences_[i].bearing_vector_0;
                    residual(i) = (pow((v0.transpose() * E * v1), 2));
            }
            return residual;
        }        
        
        
        // update points with weights
        void EssentialProblem::updateWeights(const weights_t & new_weights)
        {

            assert(new_weights.cols() == number_points_);
            for(size_t i = 0; i < number_points_; i++)   point_correspondences_[i].weight_ = new_weights(i);


            // update data matrix
            data_matrix_C_ = constructDataMatrix(point_correspondences_);
        }

        Matrix34 EssentialProblem::initializeSolver8pts(void)
        {
          // we suppose that data_matrix_C_ is already updated!!
          // for simplicity, we only implement the 8pt algorithm
          // Compute E AND save precon
          Matrix3 E = initialize8pts(data_matrix_C_, Matrix_precon_, Preconditioner::Any);
          
          // point_correspondences_ are already updated!!
          Matrix3 R_s;
          Vector3 t_s;
          computeRtfromE(point_correspondences_, E, R_s, t_s);
          Matrix34 Rt;
          Rt.block<3, 3>(0, 0) = R_s;
          Rt.block<3, 1>(0, 3) = t_s;
          return (Rt);
        }

       

        // Compute pseudo Jacobi preconditioner
        Matrix3 EssentialProblem::computePseudoJacobiPrecon(void)
        {
                // If you use the 8 points algorithm, please
                // select the use_precon option

                Eigen::JacobiSVD<Matrix9> svd(data_matrix_C_, Eigen::ComputeFullU | Eigen::ComputeFullV);  //faster

                Vector9 svd_diag = svd.singularValues();
                return (Matrix3::Identity() / svd_diag(0));
        }

        Matrix3 EssentialProblem::invMtPrecon(Matrix3 & R_init)
        {
                Matrix3 Mt = createMatrixT(data_matrix_C_, R_init);
                Eigen::JacobiSVD<Matrix3> svd(Mt, Eigen::ComputeFullU | Eigen::ComputeFullV);
                // return (Mt);
                Vector3 svd_diag = svd.singularValues();
                return (Matrix3::Identity() / svd_diag(0));
        }


      

        
        // apply precon to the full X
        Matrix34 EssentialProblem::precondition(const Matrix34& X, const Matrix34 & Xdot) const
        {
             if (use_precon_ == Preconditioner::None)
                        return Xdot;
             else       return tangent_space_projection(X,  Matrix_precon_ * Xdot);
        }
        
        
     


        double EssentialProblem::evaluate_objective(const Matrix34 &Y) const {
            // TODO: save this for latter (gradient & hessian)
            Vector3 t = Y.block<3, 1>(0, 3);
            Matrix3 R = Y.block<3, 3>(0, 0);
            Matrix3 Mt = createMatrixT(data_matrix_C_, R);
            
            Matrix9 Mr = createMatrixR(data_matrix_C_, t);
            

            return (0.5 * (t.transpose() * Mt * t).trace());

        }

         double EssentialProblem::evaluate_objective(const Matrix34 &Y, ProblemCachedMatrices & problem_matrices) const {

            Vector3 t = Y.block<3, 1>(0, 3);
            Matrix3 R = Y.block<3, 3>(0, 0);
            problem_matrices.Mt = createMatrixT(data_matrix_C_, R);
            problem_matrices.Mr = createMatrixR(data_matrix_C_, t);

            return (0.5 * (t.transpose() * problem_matrices.Mt * t).trace());

        }


            // TODO: Use the fcn below
            Matrix34 EssentialProblem::Euclidean_gradient(const Matrix34 &Y) const
            {
                Vector3 t = Y.block<3, 1>(0, 3);
                Matrix3 R = Y.block<3, 3>(0, 0);
                Matrix3 Mt = createMatrixT(data_matrix_C_, R);
                Matrix9 Mr = createMatrixR(data_matrix_C_, t);
                Matrix34 G;
                G.setZero();
                Vector9 mr = Mr * vec(R);
                G.block<3,3>(0,0)= Eigen::Map<Matrix3> (mr.data(), 3, 3);
                G.block<3,1>(0,3) = Mt * t;
                return G;
            }

             Matrix34 EssentialProblem::Euclidean_gradient(const Matrix34 &Y, const ProblemCachedMatrices & problem_matrices) const
            {
                Vector3 t = Y.block<3, 1>(0, 3);
                Matrix3 R = Y.block<3, 3>(0, 0);
                Matrix34 G;
                G.setZero();
                Vector9 mr = problem_matrices.Mr * vec(R);
                G.block<3,3>(0,0)= Eigen::Map<Matrix3> (mr.data(), 3, 3);
                G.block<3,1>(0,3) = problem_matrices.Mt * t;
                return G;
            }


            Matrix34 EssentialProblem::Riemannian_gradient(const Matrix34 &Y, const Matrix34 &nablaF_Y) const
             {
              return tangent_space_projection(Y, nablaF_Y);
             }

             Matrix34 EssentialProblem::Riemannian_gradient(const Matrix34 &Y) const {
              return tangent_space_projection(Y, Euclidean_gradient(Y));
            }

            Matrix34 EssentialProblem::Riemannian_gradient(const Matrix34 &Y, const ProblemCachedMatrices & problem_matrices) const {
              return tangent_space_projection(Y, Euclidean_gradient(Y, problem_matrices));
            }


       /** Given a matrix Y in the domain D of the SE-Sync optimization problem, and
           * a tangent vector dotY in T_D(Y), the tangent space of the domain of the
           * optimization problem at Y, this function computes and returns Hess
           * F(Y)[dotY], the action of the Riemannian Hessian on dotX */
           // TODO: same with gradient
               Matrix34 EssentialProblem::Riemannian_Hessian_vector_product(const Matrix34 &Y,
                                                   const Matrix34 &nablaF_Y,
                                                   const Matrix34 &dotY) const
                                                   {
                                                   // Euclidean Hessian-vector product
                                                    Matrix34 HessRiemannian;

                                                    Vector3 t = Y.block<3, 1>(0, 3);
                                                    Matrix3 R = Y.block<3, 3>(0, 0);


                                                    Matrix3 Mt = createMatrixT(data_matrix_C_, R);
                                                    Matrix9 Mr = createMatrixR(data_matrix_C_, t);
                                                    Vector3 Vt = dotY.block<3, 1>(0, 3);
                                                    Matrix3 VR = dotY.block<3, 3>(0, 0);

                                                    Matrix39 mixed_terms_der = constructMixedDifTtr(B1_,
                                                                 B2_, B3_, data_matrix_C_ , t, R);

                                                    // Compute Euclidean Hessian
                                                    Matrix9By12 coeff_R;
                                                    coeff_R.block<9, 9>(0, 0) = Mr;
                                                    coeff_R.block<9, 3>(0, 9) = 2 * mixed_terms_der.transpose();

                                                    Vector12 Vx;
                                                    Vx.block<9, 1>(0, 0) = vec(VR);
                                                    Vx.block<3, 1>(9, 0) = Vt;

                                                    Vector9 mr = coeff_R * Vx;
                                                    Matrix3 HessR = Eigen::Map<Matrix3> (mr.data(), 3, 3);

                                                    Matrix3By12 coeff_T;
                                                    coeff_T.block<3, 3>(0, 9) = Mt;
                                                    coeff_T.block<3, 9>(0, 0) = 2 * mixed_terms_der;

                                                    Vector3 HessT = coeff_T * Vx;

                                                    // recover NabldaF(Y)
                                                    Vector3 nabla_Yt = nablaF_Y.block<3,1>(0,3);
                                                    Matrix3 nabla_YR = nablaF_Y.block<3,3>(0,0);

                                                    // clean
                                                    HessRiemannian.setZero();

                                                    // Riemannain Hessian for t (sphere)
                                                    HessRiemannian.block<3,1>(0, 3) = domain_.ProjSphere(t, HessT) - (t.dot(nabla_Yt)) * Vt;

                                                    // Riemannain Hessian for R (rotation)
                                                    HessRiemannian.block<3,3>(0, 0) = domain_.ProjRotation(R, HessR - domain_.SymProduct(VR, R, nabla_YR));
                                                    // Riemannian hessian
                                                    return HessRiemannian;

                                                   }

     Matrix34 EssentialProblem::Riemannian_Hessian_vector_product(const Matrix34 &Y,
                                                   const ProblemCachedMatrices & problem_matrices,
                                                   const Matrix34 &dotY) const
                                                   {
                                                   // Euclidean Hessian-vector product
                                                    Matrix34 HessRiemannian;

                                                    Vector3 t = Y.block<3, 1>(0, 3);
                                                    Matrix3 R = Y.block<3, 3>(0, 0);


                                                    Vector3 Vt = dotY.block<3, 1>(0, 3);
                                                    Matrix3 VR = dotY.block<3, 3>(0, 0);

                                                    Matrix39 mixed_terms_der = constructMixedDifTtr(B1_,
                                                                 B2_, B3_, data_matrix_C_ , t, R);

                                                    // Compute Euclidean Hessian
                                                    Matrix9By12 coeff_R;
                                                    coeff_R.block<9, 9>(0, 0) = problem_matrices.Mr;
                                                    coeff_R.block<9, 3>(0, 9) = 2 * mixed_terms_der.transpose();

                                                    Vector12 Vx;
                                                    Vx.block<9, 1>(0, 0) = vec(VR);
                                                    Vx.block<3, 1>(9, 0) = Vt;

                                                    Vector9 mr = coeff_R * Vx;
                                                    Matrix3 HessR = Eigen::Map<Matrix3> (mr.data(), 3, 3);

                                                    Matrix3By12 coeff_T;
                                                    coeff_T.block<3, 3>(0, 9) = problem_matrices.Mt;
                                                    coeff_T.block<3, 9>(0, 0) = 2 * mixed_terms_der;

                                                    Vector3 HessT = coeff_T * Vx;

                                                    // recover NabldaF(Y)
                                                    Vector3 nabla_Yt = (problem_matrices.NablaF_Y).block<3,1>(0,3);
                                                    Matrix3 nabla_YR = (problem_matrices.NablaF_Y).block<3,3>(0,0);

                                                    // clean for sanity
                                                    HessRiemannian.setZero();

                                                    // Riemannain Hessian for t (sphere)
                                                    HessRiemannian.block<3,1>(0, 3) = domain_.ProjSphere(t, HessT) - (t.dot(nabla_Yt)) * Vt;

                                                    // Riemannain Hessian for R (rotation)
                                                    HessRiemannian.block<3,3>(0, 0) = domain_.ProjRotation(R, HessR - domain_.SymProduct(VR, R, nabla_YR));
                                                    // Riemannian hessian
                                                    return HessRiemannian;
                           }


          Matrix34 EssentialProblem::Riemannian_Hessian_vector_product(const Matrix34 &Y,
                                                    const Matrix34 &dotY) const
                                                    { return Riemannian_Hessian_vector_product(Y, Euclidean_gradient(Y), dotY);}



          Matrix34 EssentialProblem::tangent_space_projection(const Matrix34 &Y,
                                               const Matrix34 &dotY) const { return domain_.Proj(Y, dotY); }


            Matrix34 EssentialProblem::retract(const Matrix34 &Y, const Matrix34 &dotY) const
            {
                return domain_.retract(Y, dotY);

            }

          Matrix34 EssentialProblem::random_sample() const
          {
            return domain_.random_sample();
          }

} // end of Essential namespace
