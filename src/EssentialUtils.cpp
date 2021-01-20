#include "EssentialUtils.h"

#include <set>
#include <stdlib.h>  /* srand, rand */

namespace Essential{

        Eigen::MatrixXf symmetrize(const Eigen::MatrixXf & matrix_data)
        {  return 0.5 * (matrix_data + matrix_data.transpose());   }


        Vector9 vec(Matrix3 & M)
        {    return Eigen::Map<Vector9> (M.data(), 9, 1);       }

        Matrix3 cross(const Vector3 & t)
        {
                Matrix3 t_cross;
                t_cross.setZero();
                t_cross(0, 1) = -t(2); t_cross(0, 2) = t(1);
                t_cross(1, 0) = t(2); t_cross(1, 2) = -t(0);
                t_cross(2, 0) = -t(1); t_cross(2, 1) = t(0);
                return t_cross;
        }

        Eigen::MatrixXf skew(const Eigen::MatrixXf & M)
        {
            return (0.5 * (M - M.transpose()));
        }


        Matrix39 createRhat(const Matrix3& R)
        {
            Matrix39 R_hat;
            R_hat.setZero();

            for (int i = 0; i < 3; i++)
            {
                Vector3 column_i = R.col(i);
                Matrix3 temp_cross = cross(-column_i);
                R_hat.block<3, 3>(0, i*3) = temp_cross.transpose();  

            }
            return R_hat;
        }

        Matrix9 createThat(const Vector3 & t)
        {
                Matrix9 that; that.setZero();
                Matrix3 temp_cross;
                temp_cross.setZero();
                temp_cross = cross(t);

                for (int i = 0; i < 3; i++) that.block<3, 3>(i*3, i*3) = temp_cross;
                return (that);
        }

        Matrix3 createMatrixT(const Matrix9 & C, const Matrix3 & R)
        {
                const Matrix39 R_hat = createRhat(R) ;
                const Matrix39 temp = R_hat * C;
                const Matrix3 temp_final = temp * (R_hat.transpose());

                return (0.5 * (temp_final + temp_final.transpose()));
        }

        Matrix9 createMatrixR(const Matrix9 & C, const Vector3 & t)
        {

                Matrix9 that = createThat(t);
                Matrix9 temp = (that.transpose() * (C * that));
                return 0.5 * (temp + temp.transpose());

        }

        Matrix9 constructDataMatrix(const bearing_vectors_t & bearing_vectors)
        {

                Matrix9 C;
                Matrix9 temp;
                // clean output matrix
                C.setZero();

                for (int i = 0; i < bearing_vectors.size(); i++)
                {
                        // clean data
                        temp.setZero();
                        const Vector3 v1 = bearing_vectors[i].bearing_vector_1;
                        const Vector3 v0 = bearing_vectors[i].bearing_vector_0;
                        const double weight = bearing_vectors[i].weight_;
                        for (int j = 0; j < 3; j++)    temp.block<3, 1>(j*3, 1) = v1[j] * v0;

                        C += weight * temp * temp.transpose();
                }
               
                return 0.5 * (C + C.transpose());
        }

        Matrix9 constructDataMatrix(const bearing_vectors_t & bearing_vectors, const weights_t & new_weights)
        {

          Matrix9 C;
          Matrix9 temp;
          // clean output matrix
          C.setZero();

          for (int i = 0; i < bearing_vectors.size(); i++)
          {
                  // clean data
                  temp.setZero();
                  const Vector3 v1 = bearing_vectors[i].bearing_vector_1;
                  const Vector3 v0 = bearing_vectors[i].bearing_vector_0;
                  const double weight = new_weights[i];
                  for (int j = 0; j < 3; j++)    temp.block<3, 1>(j*3, 1) = v1[j] * v0;
                  C += weight * temp * temp.transpose();
          }
          
          return 0.5 * (C + C.transpose());
        }



        Matrix39 constructMixedDifTtr(const Matrix9 & B1, const Matrix9 & B2,
                const Matrix9 & B3,const Matrix9 & C , const Vector3 & t, Matrix3 & R)
        {
            Matrix39 mixed_t;
            Matrix9 temp_coeff;
            Vector9 term_ii;


            Matrix9 that = createThat(t);
            
            // NOTE: avoid the loop and define it directly
            // clean everything
            term_ii.setZero();
            mixed_t.setZero();
            temp_coeff.setZero();

            // first row
            temp_coeff = B1 * C * that;
            term_ii = ( 0.5 * (temp_coeff + temp_coeff.transpose())) * vec(R);
            mixed_t.block<1, 9>(0, 0) = term_ii.transpose();

            // second row
            term_ii.setZero();
            temp_coeff.setZero();
            temp_coeff = B2 * C * that;
            term_ii = ( 0.5 * (temp_coeff + temp_coeff.transpose())) * vec(R);
            mixed_t.block<1, 9>(1, 0) = term_ii.transpose();

            // third row
            term_ii.setZero();
            temp_coeff.setZero();
            temp_coeff = B3 * C * that;
            term_ii = ( 0.5 * (temp_coeff + temp_coeff.transpose())) * vec(R);
            mixed_t.block<1, 9>(2, 0) = term_ii.transpose();

           return mixed_t;
        }

        void constructBMatricesForMixedDiftR(Matrix9& B1, Matrix9 & B2, Matrix9 & B3)
        {
            // useful vars
            size_t r1 = 0, r4 =1, r7 = 2, r2 = 3, r5 = 4, r8 = 5, r3 = 6, r6 = 7, r9 = 8;

            // create first matrix
            B1.setZero();
            B1(1, r7) = -1; B1(2, r4) = 1; B1(4, r8) = -1; B1(5, r5) = 1; B1(7, r9) = -1; B1(8, r6) = 1;
            B1.transposeInPlace();

            // create second matrix
            B2.setZero();
            B2(0, r7) = 1; B2(2, r1) = -1; B2(3, r8) = 1; B2(5, r2) = -1; B2(6, r9) = 1; B2(8, r3) = -1;
            B2.transposeInPlace();

            // create third matrix
            B3.setZero();
            B3(0, r4) = -1; B3(1, r1) = 1; B3(3, r5) = -1; B3(4, r2) = 1; B3(6, r6) = -1; B3(7, r3) = 1;
            B3.transposeInPlace();

        }


        /*
        Function adapted from opengv
        see opengv::point_t opengv::triangulation::triangulate2
        */
        Vector3 triangulatePoint(const Matrix3 & R, const Vector3 & t, const Vector3 & f0, const Vector3& f1)
        {
          // create poses
          Matrix34 P1 = Matrix34::Zero();
          Matrix34 P2 = Matrix34::Zero();

          P1.block<3,3>(0, 0) = Matrix3::Identity();
          P2.block<3,3>(0, 0) = R.transpose();
          P2.block<3, 1>(0, 3) = R.transpose() * t;

          Matrix4 A = Matrix4::Zero();
          A.row(0) = f0[0] * P1.row(2) - f0[2] * P1.row(0);
          A.row(1) = f0[1] * P1.row(2) - f0[2] * P1.row(1);
          A.row(2) = f1[0] * P2.row(2) - f1[2] * P2.row(0);
          A.row(3) = f1[1] * P2.row(2) - f1[2] * P2.row(1);
          Eigen::JacobiSVD< Eigen::MatrixXd > SVD_A(A, Eigen::ComputeFullV );


          Vector3 world_point;
          world_point[0] = SVD_A.matrixV()(0,3);
          world_point[1] = SVD_A.matrixV()(1,3);
          world_point[2] = SVD_A.matrixV()(2,3);
          world_point = world_point / SVD_A.matrixV()(3,3);

          return world_point;

        }

        Matrix3 projectToRotation(const Matrix3 & R_hat)
        {
                /*
                - Find the closest matrix under Frobenius norm
                - Use the SVD decomposition of R_hat
                1. R_hat = US*V'
                2. If det(U * V') > 0: R = U*V';
                else                 : R = U * diag(1, 1, -1) * V';

                */
                const Eigen::JacobiSVD<Matrix3> svd(R_hat, Eigen::ComputeFullU | Eigen::ComputeFullV);

                Matrix3 U = svd.matrixU();
                Matrix3 Vt = (svd.matrixV()).transpose();

                if ((U*Vt).determinant() > 0)  return (U*Vt);
                else
                {
                        // Negative determinant!
                        const Eigen::DiagonalMatrix<double, 3> d(1, 1, -1);
                        return (U*d*Vt);
                }

        }


        Matrix3 projectToEssentialManifold(const Matrix3 & E_hat)
        {
            const Eigen::JacobiSVD<Matrix3> svd(E_hat, Eigen::ComputeFullU | Eigen::ComputeFullV);
            
            const Eigen::DiagonalMatrix<double, 3> d(1, 1, 0);  // imply: ||E||_f ^2 = 2

            // check determinant
            Matrix3 Ur = svd.matrixU();
            Matrix3 Vrt = svd.matrixV().transpose();

            return (Ur * d * Vrt);
        }

        // Basic
        Matrix3 initialize8pts(const bearing_vectors_t & points)
        {
            Matrix9 C;
            Matrix3 precon;

            C = constructDataMatrix(points);
            Matrix3 E_initial = initialize8pts(C, precon, Preconditioner::None);
            return E_initial;
        }


        // 8 points:
        Matrix3 initialize8pts(const bearing_vectors_t & points, Matrix9 & C, Matrix3 & matrix_precon, Preconditioner use_precon)

        {
            C = constructDataMatrix(points);
            Matrix3 E_initial = initialize8pts(C, matrix_precon, use_precon);
            return E_initial;
        }


        Matrix3 initialize8pts(const Matrix9 & C, Matrix3 & matrix_precon, Preconditioner use_precon)
        {
            Vector9 vecE;
            vecE.setZero();
            Matrix3 E_initial;
      
                Eigen::JacobiSVD<Matrix9> svd(C, Eigen::ComputeFullU | Eigen::ComputeFullV);  //faster
                // Eigen::BDCSVD<Matrix9> svd(C, Eigen::ComputeFullU | Eigen::ComputeFullV);
                Vector9 sigmas = svd.singularValues();


                if (use_precon == Preconditioner::N)
                        matrix_precon = (Matrix3::Identity() / (sigmas.sum()) );  // = N

                else if (use_precon == Preconditioner::Dominant_eigenvalues)
                        matrix_precon = (Matrix3::Identity() / (sigmas(0) + sigmas(1) + sigmas(2)) );
                
                else
                        matrix_precon = (Matrix3::Identity() / sigmas(0));
                

                vecE = svd.matrixV().col(8);
            

            E_initial = Eigen::Map<Matrix3>(vecE.data(), 3, 3);

            return (projectToEssentialManifold(E_initial));
        }


        /* Minimal solvers */
        // 8-pt
        Matrix3 minimalinitialize8ptsMinimal(const bearing_vectors_t & points)
        {
            // Create adapter for opengv
            std::vector<opengv::bearingVector_t, Eigen::aligned_allocator<opengv::bearingVector_t> > bearingVectors1, bearingVectors2;

            double N = points.size();
            // use only eight correspondences
            std::vector<int> set_idx =  selectKPoints(N, 8);
            
                for (int i = 0; i < set_idx.size(); i++)
                    {
                        opengv::bearingVector_t p1 = points[set_idx[i]].bearing_vector_0.normalized();
                        opengv::bearingVector_t p2 = points[set_idx[i]].bearing_vector_1.normalized();

                        bearingVectors1.push_back(p1);
                        bearingVectors2.push_back(p2);
                    }

            opengv::relative_pose::CentralRelativeAdapter adapter( bearingVectors1, bearingVectors2, Eigen::Matrix3d::Identity() );
            
            auto E8 = opengv::relative_pose::eightpt( adapter );
            
            // Call actual 8 points
            return (projectToEssentialManifold(E8));
        }
        
        
        
        
        Matrix3 minimalinitialize8ptsAll(const bearing_vectors_t & points)
        {
            // Create adapter for opengv
            std::vector<opengv::bearingVector_t, Eigen::aligned_allocator<opengv::bearingVector_t> > bearingVectors1, bearingVectors2;

            double N = points.size();
            // Here we use all the correspondences, i.e., DLT
                for (int i = 0; i < N; i++)
                    {
                        opengv::bearingVector_t p1 = points[i].bearing_vector_0.normalized();
                        opengv::bearingVector_t p2 = points[i].bearing_vector_1.normalized();

                        bearingVectors1.push_back(p1);
                        bearingVectors2.push_back(p2);
                    }

            opengv::relative_pose::CentralRelativeAdapter adapter( bearingVectors1, bearingVectors2, Eigen::Matrix3d::Identity() );
            
            auto E8 = opengv::relative_pose::eightpt( adapter );
            
            // Call actual 8 points
            return (projectToEssentialManifold(E8));
        }
        
        
        Matrix3 minimalinitialize8pts(const bearing_vectors_t & points, bool use_all_corr)
        {
                if (use_all_corr)
                     return minimalinitialize8ptsAll(points);   
        
                else 
                     return minimalinitialize8ptsMinimal(points); 
        }

        
        
        /* Select rnd K points without repetition in 1, ..., max_n */
        std::vector<int> selectKPoints(int max_n, int k)
        {
                std::vector<int> out; 
                std::set<int> set_n;
        
                int idx = 0; 
                bool new_idx = true;
                
                for (int i = 0; i < k; i++)
                {
                        do{
                                idx = std::round( (max_n - 1) * ( (double) rand() ) / ( (double) RAND_MAX ) );
                                
                                // if element is not in the set
                                if (set_n.count(idx) == 0)
                                {
                                        out.push_back(idx); 
                                        set_n.insert(idx); 
                                        new_idx = true;  // break loop
                                }
                                else new_idx = false;
                
                        }while(!new_idx);
                
                }
                
                return out;
        
        }


        // 5-pt algorithm
        Matrix3 initialize5ptsMinimal(const bearing_vectors_t & points)
        {
                // Create adapter for opengv
                std::vector<opengv::bearingVector_t, Eigen::aligned_allocator<opengv::bearingVector_t> > bearingVectors1, bearingVectors2;

                double N = points.size();
                // use only five points 
                std::vector<int> set_idx =  selectKPoints(N, 5);
                
                for (int i = 0; i < set_idx.size(); i++)
                    {
                        opengv::bearingVector_t p1 = points[set_idx[i]].bearing_vector_0.normalized();
                        opengv::bearingVector_t p2 = points[set_idx[i]].bearing_vector_1.normalized();

                        bearingVectors1.push_back(p1);
                        bearingVectors2.push_back(p2);
                    }

                opengv::relative_pose::CentralRelativeAdapter adapter( bearingVectors1, bearingVectors2, Eigen::Matrix3d::Identity() );

                // Call actual 5 points
                Matrix3 E_5pt;
                opengv::essentials_t fivept_nister_essentials = opengv::relative_pose::fivept_nister( adapter );
                double sum_res = 0.0, min_res = 100000.0;
                // for each solution
                size_t N_red=4;
                
                for( size_t i = 0; i < fivept_nister_essentials.size(); i++ )
                {
                  // compute the residuals and sum them up
                  sum_res = 0.0;
                  for (size_t j = 0; j < N_red; j++)
                  {
                    Vector3 v1 = points[j].bearing_vector_1;
                    Vector3 v0 = points[j].bearing_vector_0;
                    sum_res += (pow((v0.transpose() * fivept_nister_essentials.at(i) * v1), 2));
                  }
                  
                  if (sum_res < min_res)
                  {
                    min_res = sum_res;  // update min. residual
                    E_5pt = fivept_nister_essentials.at(i);  // save current solution
                  }
                }

                return (projectToEssentialManifold(E_5pt));
        }
        
        
        
        
        
        Matrix3 initialize5ptsAll(const bearing_vectors_t & points)
        {
               // Create adapter for opengv
                std::vector<opengv::bearingVector_t, Eigen::aligned_allocator<opengv::bearingVector_t> > bearingVectors1, bearingVectors2;

                double N = points.size();
                for (int i = 0; i < N; i++)
                    {
                        opengv::bearingVector_t p1 = points[i].bearing_vector_0.normalized();
                        opengv::bearingVector_t p2 = points[i].bearing_vector_1.normalized();

                        bearingVectors1.push_back(p1);
                        bearingVectors2.push_back(p2);
                    }

                opengv::relative_pose::CentralRelativeAdapter adapter( bearingVectors1, bearingVectors2, Eigen::Matrix3d::Identity() );

                // Call actual 5 points
                Matrix3 E_5pt;
                opengv::essentials_t fivept_nister_essentials = opengv::relative_pose::fivept_nister( adapter );
                double sum_res = 0.0, min_res = 100000.0;
                // for each solution
                size_t N_red=4;
                
                for( size_t i = 0; i < fivept_nister_essentials.size(); i++ )
                {
                  // compute the residuals and sum them up
                  sum_res = 0.0;
                  for (size_t j = 0; j < N_red; j++)
                  {
                    Vector3 v1 = points[j].bearing_vector_1;
                    Vector3 v0 = points[j].bearing_vector_0;
                    sum_res += (pow((v0.transpose() * fivept_nister_essentials.at(i) * v1), 2));
                  }
                  if (sum_res < min_res)
                  {
                    min_res = sum_res;  // update min. residual
                    E_5pt = fivept_nister_essentials.at(i);  // save current solution
                  }
                  
                }

                return (projectToEssentialManifold(E_5pt));
        }
        
        
        
        
        Matrix3 initialize5pts(const bearing_vectors_t & points, bool use_all_corr)
        {
                if (use_all_corr)
                        return initialize5ptsAll(points);
                else 
                        return initialize5ptsMinimal(points);        
        }
        
        

        void computeRtfromE(const bearing_vectors_t & points, const Matrix3& E, Matrix3 & R, Vector3 & t )
        {
            // Note: E should be already in \Me

            Eigen::JacobiSVD<Matrix3> svd(E, Eigen::ComputeFullU | Eigen::ComputeFullV);

            Matrix3 W;
            W.setZero();
            W(0, 1) = -1; W(1, 0) = 1; W(2, 2) = 1;

            Matrix3 R1, R2, V, U;
            U = svd.matrixU();
            V = svd.matrixV();

            // Rotation matrix
            R1 = U * W * V.transpose();

            R2 = U * W.transpose() * V.transpose();

            // for R1, R2 in SO(3)
            if (R1.determinant() < 0)       R1 = -R1;
            if (R2.determinant() < 0)       R2 = -R2;


            // Translation vector
            Vector3 t1;
            t1 = U.col(2);


            // Cheirality

            unsigned int N = points.size();
            Eigen::MatrixXd X(9, N);
            Eigen::MatrixXd Y(6, N);
            X.setZero();
            Y.setZero();


            for (int i = 0; i < N; i++)
            {
                Vector9 temp;
                Vector3 v1 = points[i].bearing_vector_1;
                Vector3 v0 = points[i].bearing_vector_0;
                double weight_i = points[i].weight_;

                temp.setZero();
                for (int j = 0; j < 3; j++)    temp.block<3, 1>(j*3, 0) = weight_i * v1[j] * v0;
                X.col(i) = temp;


                Y.block<3, 1>(0, i) = v0 * v1.norm(); Y.block<3, 1>(3, i) = v1 * v0.norm();

            }

            Matrix3 ER1 = E * E.transpose() * R1;
            Matrix3 ER2 = E * E.transpose() * R2;

            Eigen::MatrixXd M11(N, 1);
            Eigen::MatrixXd M12(N, 1);
            M11 = X.transpose() * Eigen::Map<Eigen::Matrix<double, 9, 1> > (ER1.data(), 9, 1);
            M12 = X.transpose() * Eigen::Map<Eigen::Matrix<double, 9, 1> > (ER2.data(), 9, 1);
            double n_non_zeros_M11 = (M11.array() > 0).count();
            double n_non_zeros_M12 = (M12.array() > 0).count();


            if (n_non_zeros_M11 >= n_non_zeros_M12)      R = R1;
            else                                         R = R2;


            // Compute the right translation

            Eigen::MatrixXd M21(N, 1), M22(N, 1);
            Eigen::MatrixXd R_eye(6, 3); R_eye.block<3,3>(0,0) = - R.transpose(); R_eye.block<3,3>(3, 0) = Matrix3::Identity();


            M21 = Y.transpose() * R_eye * t1;

            double n_non_zeros_M21 = (M21.array() > 0).count();
            double n_non_zeros_M22 = (-M21.array() > 0).count();


            if (n_non_zeros_M21 >= n_non_zeros_M22)    t = t1;
            else                                       t = -t1;

        };



         Matrix3 initialize7pts(const bearing_vectors_t & points, bool use_all_corr)
         {
                if (use_all_corr)
                        return initialize7ptsAll(points);
                else 
                        return initialize7ptsMinimal(points);
         
         
         }
         
         
         Matrix3 initialize7ptsMinimal(const bearing_vectors_t & points)
         {
                // Create adapter for opengv
                std::vector<opengv::bearingVector_t, Eigen::aligned_allocator<opengv::bearingVector_t> > bearingVectors1, bearingVectors2;

                double N = points.size();
                // use only five points 
                std::vector<int> set_idx =  selectKPoints(N, 7);
                
                for (int i = 0; i < set_idx.size(); i++)
                    {
                        opengv::bearingVector_t p1 = points[set_idx[i]].bearing_vector_0.normalized();
                        opengv::bearingVector_t p2 = points[set_idx[i]].bearing_vector_1.normalized();

                        bearingVectors1.push_back(p1);
                        bearingVectors2.push_back(p2);
                    }

                opengv::relative_pose::CentralRelativeAdapter adapter( bearingVectors1, bearingVectors2, Eigen::Matrix3d::Identity() );

                // Call actual 7 points
                opengv::essentials_t seven_pt_essentials = opengv::relative_pose::sevenpt(adapter);

                // obtain the solution
                Matrix3 E_7pt;
                double sum_res = 0.0, min_res = 100000.0;
                size_t N_red = 4;
                for (size_t j = 0; j < seven_pt_essentials.size(); j++)
                {
                  for (int i = 0; i < N_red; i++)
                  {
                    Vector3 v1 = points[i].bearing_vector_1;
                    Vector3 v0 = points[i].bearing_vector_0;
                    sum_res += (pow((v0.transpose() * seven_pt_essentials.at(j) * v1), 2));
                  }
                  if (sum_res < min_res)
                  {
                    min_res = sum_res;  // update min. residual
                    E_7pt = seven_pt_essentials.at(j);  // save current solution
                  }
                  sum_res = 0.0;
                }

                return (projectToEssentialManifold(E_7pt));
         }
         
         
         
         
          Matrix3 initialize7ptsAll(const bearing_vectors_t & points)
         {
                // Create adapter for opengv
                std::vector<opengv::bearingVector_t, Eigen::aligned_allocator<opengv::bearingVector_t> > bearingVectors1, bearingVectors2;

                double N = points.size(); 
                
                for (int i = 0; i < N; i++)
                    {
                        opengv::bearingVector_t p1 = points[i].bearing_vector_0.normalized();
                        opengv::bearingVector_t p2 = points[i].bearing_vector_1.normalized();

                        bearingVectors1.push_back(p1);
                        bearingVectors2.push_back(p2);
                    }

                opengv::relative_pose::CentralRelativeAdapter adapter( bearingVectors1, bearingVectors2, Eigen::Matrix3d::Identity() );

                // Call actual 7 points
                opengv::essentials_t seven_pt_essentials = opengv::relative_pose::sevenpt(adapter);

                // obtain the solution
                Matrix3 E_7pt;
                double sum_res = 0.0, min_res = 100000.0;
                size_t N_red = 4;
                for (size_t j = 0; j < seven_pt_essentials.size(); j++)
                {
                  for (int i = 0; i < N_red; i++)
                  {
                    Vector3 v1 = points[i].bearing_vector_1;
                    Vector3 v0 = points[i].bearing_vector_0;
                    sum_res += (pow((v0.transpose() * seven_pt_essentials.at(j) * v1), 2));
                  }
                  if (sum_res < min_res)
                  {
                    min_res = sum_res;  // update min. residual
                    E_7pt = seven_pt_essentials.at(j);  // save current solution
                  }
                  sum_res = 0.0;
                }

                return (projectToEssentialManifold(E_7pt));
         }
         
         
         
         


        bearing_vectors_t convertPointsVector32PointsOpt(std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> & bearingVectors0,
                                                    std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> & bearingVectors1 )
        {
            bearing_vectors_t points_opt;

            double N = bearingVectors1.size();

            for (int i = 0; i < N; i++)
            {
                Vector3 p0 = bearingVectors0[i] / bearingVectors0[i].norm();
                Vector3 p1 = bearingVectors1[i] / bearingVectors1[i].norm();
                double weighti = 1.0;

                CorrespondingFeatures bearing(p0, p1, weighti);
                points_opt.push_back(bearing);
            }
            return points_opt;

        }

        bearing_vectors_t convertPointsVector32PointsOpt(std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> & bearingVectors0,
                                                    std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> & bearingVectors1,
                                                    std::vector<double> & weights_observations)
        {
            bearing_vectors_t points_opt;
            // TODO Assert dim P1, P2, and W
            double N = bearingVectors1.size();

            for (int i = 0; i < N; i++)
            {
                Vector3 p0 = bearingVectors0[i] / bearingVectors0[i].norm();
                Vector3 p1 = bearingVectors1[i] / bearingVectors1[i].norm();
                double weighti = weights_observations[i];

                CorrespondingFeatures bearing(p0, p1, weighti);
                points_opt.push_back(bearing);
            }
            return points_opt;

        }


        bearing_vectors_t convertPointsVector32PointsOpt(Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors0,
                                                         Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors1,
                                                         Eigen::Matrix<double, 1, Eigen::Dynamic> & weights_observations)
        {
            bearing_vectors_t points_opt;
            // TODO Assert dim P1, P2, and W
            double N = bearingVectors1.cols();

            for (int i = 0; i < N; i++)
            {
                Vector3 p0 = bearingVectors0.col(i) / bearingVectors0.col(i).norm();
                Vector3 p1 = bearingVectors1.col(i) / bearingVectors1.col(i).norm();
                double weighti = weights_observations(i);

                CorrespondingFeatures bearing(p0, p1, weighti);
                points_opt.push_back(bearing);
            }
            return points_opt;

        }
        


        Matrix3 computeEfromRt(const Matrix3 & R, const Vector3 & t)
        {
            Matrix3 t_skew = cross(t);
            Matrix3 E = cross(t) * R;
            return E / (E.norm()) * SQRT2;
        }

        Matrix3 computeEfromRt(const Matrix34& Rt)
        {
            Vector3 t = Rt.block<3, 1>(0, 3);
            Matrix3 R = Rt.block<3, 3>(0, 0);
            return (computeEfromRt(R, t));
        }


                // compute distance (in geodesic sense) between two rotations (See paper) [degrees]
                double distR(const Matrix3 & R_gt, const Matrix3 & R_est)
                {
                    return (R180PI * acos(((R_est.transpose() * R_gt).trace() - 1) * 0.5)  );
                }
                // Compute distance (in direction) between two (unitary) translation vectors [degrees]
                double distT(const Vector3 & T_gt, const Vector3 & T_est)
                {
                    Vector3 tgt = T_gt;
                    tgt = tgt / T_gt(2); // in order to fix the sense of the direction
                    tgt.normalize();

                    Vector3 test = T_est;
                    test = test / T_est(2);
                    test.normalize();
                    return (R180PI * acos(tgt.dot(test)));
                }


                // Compute similarity between two essential matrices [degrees]
                double distE(const Matrix3 & E_gt, const Matrix3 & E_est)
                {
                    return (R180PI * acos(  (E_gt.transpose() * E_est).trace() / ( E_gt.norm() * E_est.norm() )  )  );
                }

                // Evaluate essential matrix, rotation and translation direction (everything in [degreees])

                void evaluateRTE(const Matrix3 & R_gt, const Matrix3 & R_est,
                            const Vector3 & T_gt, const Vector3 & T_est,
                            const Matrix3 & E_gt, const Matrix3 & E_est,
                            double& dist_R, double & dist_T, double & dist_E)
                            {
                                dist_R = distR(R_gt, R_est);
                                dist_T = distT(T_gt, T_est);
                                dist_E = distE(E_gt, E_est);
                            }




} //end of essential namespace
