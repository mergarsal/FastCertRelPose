#pragma once

#include <Eigen/Core>

namespace Essential{
    
                typedef Eigen::Matrix<double, 3, 1> Vector3;
                typedef Eigen::Matrix<double, 4, 1> Vector4;
                typedef Eigen::Matrix<double, 5, 1> Vector5;
                typedef Eigen::Matrix<double, 6, 1> Vector6;
                typedef Eigen::Matrix<double, 7, 1> Vector7;
                typedef Eigen::Matrix<double, 9, 1> Vector9;
                typedef Eigen::Matrix<double, 12, 1> Vector12;
                
                typedef Eigen::Matrix<double, 3, 3> Matrix3;
                typedef Eigen::Matrix<double, 4, 4> Matrix4;
                typedef Eigen::Matrix<double, 9, 9> Matrix9;
                typedef Eigen::Matrix<double, 3, 9> Matrix39;
                typedef Eigen::Matrix<double, 9, 3> Matrix93;
                typedef Eigen::Matrix<double, 3, 4> Matrix34;
                typedef Eigen::Matrix<double, 12, 12> Matrix12;
                
                typedef Eigen::Matrix<double, 9, 12> Matrix9By12;
                typedef Eigen::Matrix<double, 3, 12> Matrix3By12;
                    
               
                
                
    enum class Preconditioner {
                  None= 0,
                  Max_eigenvalue,
                  Dominant_eigenvalues, 
                  N,
                  Any,
               };
                              
} // end of essential namespace

