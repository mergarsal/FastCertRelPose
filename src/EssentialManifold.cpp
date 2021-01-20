#include "EssentialManifold.h"


namespace Essential{


            Matrix34 EssentialManifold::project(const Matrix34 & P) const
            {
                 Matrix34 A;
                 Matrix3 R_p, R = P.block<3,3>(0, 0);
                 // Classic result
                 Eigen::JacobiSVD<Matrix3> svd(R, Eigen::ComputeFullU | Eigen::ComputeFullV);

                double detU = svd.matrixU().determinant();
                double detV = svd.matrixV().determinant();

                if (detU * detV > 0)
                    R_p =  svd.matrixU() * svd.matrixV().transpose();
                else
                {
                    Matrix3 Uprime = svd.matrixU();
                    Uprime.col(Uprime.cols() - 1) *= -1;
                    R_p = Uprime * svd.matrixV().transpose();
                }

                // For the translation
                Vector3 t_p = P.block<3,1>(0, 3);
                t_p.normalize();

                // fill the output matrix
                A.setZero();
                A.block<3,3>(0, 0) = R_p;
                A.block<3,1>(0, 3) = t_p;

                return A;
            }

            Vector3 EssentialManifold::ProjSphere(const Vector3 &t, const Vector3 &Vt) const
            {  return Vt - t.dot(Vt) * t; }

            Matrix3 EssentialManifold::ProjRotation(const Matrix3& R, const Matrix3 & VR) const
            {
                    Matrix3 RtVr = R.transpose() * VR;
                    return R * 0.5 * (RtVr - RtVr.transpose());
            }

            Matrix34 EssentialManifold::retract(const Matrix34 &Y, const Matrix34 &V) const {

              // We use projection-based retraction, as described in "Projection-Like
              // Retractions on Matrix Manifolds" by Absil and Malick
              return project(Y + V);
        }
       

        Matrix34 EssentialManifold::random_sample(
            const std::default_random_engine::result_type &seed) const {
          // Generate a matrix of the appropriate dimension by sampling its elements
          // from the standard Gaussian
          std::default_random_engine generator(seed);
          std::normal_distribution<double> g;

          Matrix34 R;
          R.setZero();

          for (size_t r = 0; r < 3; ++r)
            for (size_t c = 0; c < 4; ++c)
              R(r, c) = g(generator);
          return project(R);
        }

         Matrix3 EssentialManifold::SymProduct(const Matrix3 & Ydot, const Matrix3 Y, const Matrix3 & nabla_Y) const
         {
                Matrix3 P = Y * nabla_Y;
                return Ydot.transpose() * 0.5 * (P + P.transpose());
         }

         // TODO: we can just call here to ProjSphere & ProjRotation
         Matrix34 EssentialManifold::Proj(const Matrix34 &Y, const Matrix34 &V) const{
            Matrix34 V_tan;
            V_tan.setZero();

            Matrix3 R, Vr;
            Vector3 t, Vt;

            // for the sphere
            Vt = V.block<3, 1>(0, 3);
            t = Y.block<3, 1>(0, 3);
            V_tan.block<3, 1>(0, 3) = Vt - t.dot(Vt) * t;
            // V_tan.block<3, 1>(0, 3) = ProjSphere(t, Vt);

            // for the rotation
            Vr = V.block<3,3>(0, 0);
            R = Y.block<3,3>(0, 0);
            Matrix3 RtVr = R.transpose() * Vr;
            V_tan.block<3,3>(0, 0) = R * 0.5 * (RtVr - RtVr.transpose());
            // V_tan.block<3,3>(0, 0) = ProjRotation(R, Vr);
            return V_tan;
         }


} // end of essential namespace
