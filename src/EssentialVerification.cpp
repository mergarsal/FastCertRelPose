#include "EssentialVerification.h"
#include <Eigen/Eigenvalues>



namespace Essential{

    Vector6 EssentialVerification::computeLagMultGeneral(Matrix12 & Q, size_t idx_relaxation)
    {
      /* 
      This function computes the Lagrange multipliers
      associated with the relaxation *idx_relaxation* [PARAM].
      NOTE: If you only *want* or *need* one try, use the function below
      */
      
      Vector12 x_opt, b, tt_opt;
      Vector6 Lagrange_multipliers;
      Vector5 reduced_Lagrange;
      x_opt.block<9, 1>(0, 0) = Eigen::Map<Vector9>(E_.data(), 9, 1);
      x_opt.block<3, 1>(9, 0) = t_;
      // tt_opt: blck(0_{9 times 9}, eye(3)) * x_opt
      tt_opt.setZero();
      tt_opt.block<3,1>(9,0) = t_;


      Q.setZero();
      Q.block<9, 9>(0, 0) = C_;
      // 0.5 because f = 0.5 * e' * C * e
      b = 0.5 * Q * x_opt - 0.5 * x_opt.transpose() * Q * x_opt * tt_opt;

      double x1=x_opt(0); double x2=x_opt(1); double x3=x_opt(2);
      double x4=x_opt(3); double x5=x_opt(4); double x6=x_opt(5);
      double x7=x_opt(6);double x8=x_opt(7); double x9=x_opt(8);
      double x10=x_opt(9); double x11=x_opt(10);double x12=x_opt(11);

      // int L1=0;
      int L2 = 0; int L3=L2+1; int L4=L3+1; int L5=L4+1; int L6=L5+1; int L7=L6+1;

      // NOTE; this matrix (with 6 columns) has linear dependent columns
      Eigen::Matrix<double, 12, 6> A_ext;
      A_ext.setZero();


      // first equation
      A_ext(0, L2) = x1;   A_ext(0, L3) = x2*0.5;  A_ext(0, L4) = x3*0.5;

      // Second equation
      A_ext(1, L3) = x1*0.5;  A_ext(1, L5)=x2;     A_ext(1, L6)=x3*0.5;

      // 3rd equation
      A_ext(2, L4) = x1*0.5;  A_ext(2, L6)=x2*0.5;  A_ext(2, L7) = x3;

      // 4th equation
      A_ext(3, L2)=x4;     A_ext(3, L3)=x5*0.5; A_ext(3, L4) = x6*0.5;

      // 5th equation
      A_ext(4, L3)=x4*0.5; A_ext(4, L5) = x5;   A_ext(4, L6) = x6*0.5;

      // 6th equation
      A_ext(5, L4)=x4*0.5;  A_ext(5, L6) = x5*0.5;  A_ext(5, L7) = x6;

      // 7th equation
      A_ext(6,L2)=x7;       A_ext(6, L3)=x8*0.5;  A_ext(6, L4)=x9*0.5;

      // 8th equation
      A_ext(7, L3)=x7*0.5;  A_ext(7, L5)=x8; A_ext(7, L6)=x9*0.5;

      // 9th equation
      A_ext(8, L4)=x7*0.5;  A_ext(8, L6)=x8*0.5;  A_ext(8, L7)=x9;

      // 10th equation
      A_ext(9, L3)=x11*0.5;  A_ext(9, L4)=x12*0.5;  A_ext(9, L5)=-x10;
      A_ext(9, L7)=-x10;

      // 11th equation
      A_ext(10, L3)=x10*0.5;  A_ext(10, L6)=x12*0.5; A_ext(10, L2)=-x11;
      A_ext(10, L7)=-x11;

      // 12th equation
      A_ext(11, L4)=x10*0.5;  A_ext(11, L6)=x11*0.5; A_ext(11, L2)=-x12;
      A_ext(11, L5)=-x12;

      // select relaxation
      Eigen::Matrix<double, 12, 5> A;
      A.setZero();
      // fill the matrix
      // NOTE: idx_relaxation \in {0, 5}
      // 1. check extreme cases
      if (idx_relaxation == 0)
      {
        A = A_ext.rightCols(5);
      }
      else if (idx_relaxation == 5)
      {
        A = A_ext.leftCols(5);
      }
      else
      {
        A.leftCols(idx_relaxation) = A_ext.leftCols(idx_relaxation);
        A.rightCols(5 - idx_relaxation) = A_ext.rightCols(5 - idx_relaxation);
      }

      // solving
      reduced_Lagrange = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);

      // std::cout << "Singular values lagrange:\n" << A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).singularValues() << std::endl; 

      
      // fill the vector with multipliers
      Lagrange_multipliers.setZero();
      Lagrange_multipliers(0) = 0.5 * x_opt.transpose() * Q * x_opt;
      Lagrange_multipliers.bottomRows(5) = reduced_Lagrange;


      return  Lagrange_multipliers;

    }
    
    
    
    
    

    Vector6 EssentialVerification::computeLagMult(Matrix12 & Q)
    {
      Vector12 x_opt, b, tt_opt;
      Vector6 Lagrange_multipliers;
      Vector5 reduced_Lagrange;
      x_opt.block<9, 1>(0, 0) = Eigen::Map<Vector9>(E_.data(), 9, 1);
      x_opt.block<3, 1>(9, 0) = t_;
      // tt_opt: blck(0_{9 times 9}, eye(3)) * x_opt
      tt_opt.setZero();
      tt_opt.block<3,1>(9,0) = t_;


      Q.setZero();
      Q.block<9, 9>(0, 0) = C_;
      // 0.5 because f = 0.5 * e' * C * e
      b = 0.5 * Q * x_opt - 0.5 * x_opt.transpose() * Q * x_opt * tt_opt;

      double x1=x_opt(0); double x2=x_opt(1); double x3=x_opt(2);
      double x4=x_opt(3); double x5=x_opt(4); double x6=x_opt(5);
      double x7=x_opt(6);double x8=x_opt(7); double x9=x_opt(8);
      double x10=x_opt(9); double x11=x_opt(10);double x12=x_opt(11);

      // int L1=0;
      int L3=0; int L4=1; int L5=2; int L6=3; int L7=4;
      Eigen::Matrix<double, 12, 5> A;
      A.setZero();
      // first equation
      // A(1, L2) = x1;
      A(0, L3) = x2*0.5;  A(0, L4) = x3*0.5;

      // Second equation
      A(1, L3) = x1*0.5;  A(1, L5)=x2; A(1, L6)=x3*0.5;

      // 3rd equation
      A(2, L4) = x1*0.5;  A(2, L6)=x2*0.5; A(2, L7) = x3;

      // 4th equation
      // A(4, L2) = x4;
      A(3, L3)=x5*0.5; A(3, L4) = x6*0.5;

      // 5th equation
      A(4, L3)=x4*0.5; A(4, L5) = x5;   A(4, L6) = x6*0.5;

      // 6th equation
      A(5, L4)=x4*0.5;  A(5, L6) = x5*0.5;  A(5, L7) = x6;

      // 7th equation
      // A(7, L2)=x7;
      A(6, L3)=x8*0.5;  A(6, L4)=x9*0.5;

      // 8th equation
      A(7, L3)=x7*0.5;  A(7, L5)=x8; A(7, L6)=x9*0.5;

      // 9th equation
      A(8, L4)=x7*0.5;  A(8, L6)=x8*0.5;  A(8, L7)=x9;

      // 10th equation
      A(9, L3)=x11*0.5;  A(9, L4)=x12*0.5;  A(9, L5)=-x10;
      // A(9, L1)=x10;
      A(9, L7)=-x10;

      // 11th equation
      A(10, L3)=x10*0.5;  A(10, L6)=x12*0.5;
      //A(11, L2)=-x11;
      // A(10, L1)=x11;
      A(10, L7)=-x11;

      // 12th equation
      A(11, L4)=x10*0.5;  A(11, L6)=x11*0.5;
      // A(12, L2)=-x12;
      // A(11, L1)=x12;
      A(11, L5)=-x12;

      reduced_Lagrange = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);

      // std::cout << "Singular values lagrange:\n" << A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).singularValues() << std::endl; 
      // fill the vector with multipliers
      Lagrange_multipliers.setZero();
      Lagrange_multipliers(0) = 0.5 * x_opt.transpose() * Q * x_opt;
      Lagrange_multipliers.block<5, 1>(1, 0) = reduced_Lagrange;


      return  Lagrange_multipliers;
    }
    



    double EssentialVerification::computePenalisedMatrixMAndMinEigenvalue(const Vector6 & Lagrange_multipliers, const Matrix12 & Q)
    {
        Matrix12 M;
        M.setZero();

        Matrix12 A, As;
        A.setZero();
        As.setZero();
        unsigned int e1=0; unsigned int e4=1; unsigned int e7=2;
        unsigned int e2=3; unsigned int e5=4; unsigned int e8=5;
        unsigned int e3=6; unsigned int e6=7; unsigned int e9=8;
        unsigned int t1=9; unsigned int t2=10; unsigned int t3=11;

        // initialize
        M = Q;


        // Create the matrices

        As.setZero();
        A.setZero();
        As(t1,t1)=1; As(t2,t2)=1;As(t3,t3)=1;    A = 0.5 * (As + As.transpose());
        M -=  Lagrange_multipliers(0) * A;

        As.setZero();
        A.setZero();
        As(e1,e4)=1; As(e2,e5)=1;As(e3,e6)=1;As(t1, t2)=1;  A = 0.5 * (As + As.transpose());
        M -=  Lagrange_multipliers(1) * A;

        As.setZero();
        A.setZero();
        As(e1,e7)=1; As(e2,e8)=1;As(e3,e9)=1;As(t1, t3)=1;  A = 0.5 * (As + As.transpose());
        M -=  Lagrange_multipliers(2) * A;

        As.setZero();
        A.setZero();
        As(e4,e4)=1; As(e5,e5)=1;As(e6,e6)=1;As(t1, t1)=-1; As(t3, t3)=-1; A = 0.5 * (As + As.transpose());
        M -=  Lagrange_multipliers(3) * A;

        As.setZero();
        A.setZero();
        As(e4,e7)=1; As(e5,e8)=1;As(e6,e9)=1;As(t2, t3)=1; A = 0.5 * (As + As.transpose());
        M -=  Lagrange_multipliers(4) * A;

        As.setZero();
        A.setZero();
        As(e7,e7)=1; As(e8,e8)=1;As(e9,e9)=1;As(t1, t1)=-1; As(t2, t2)=-1;  A = 0.5 * (As + As.transpose());
        M -=  Lagrange_multipliers(5) * A;


        // Compute min eigenvalue of M
        // d=eig(M);


        double min_eigenvalue = -1000;
        Vector12 min_eigenvector;



            Vector3 eigenvalues_bottomrightM = ((M.block<3,3>(9, 9)).eigenvalues()).real();
            min_eigenvalue = eigenvalues_bottomrightM.minCoeff();  
            // In many cases, this part fails
            
            if (min_eigenvalue > tau_)
            {
                // Try with the other one
                Vector9 eigenvalues_topleft = ((M.block<9,9>(0, 0)).eigenvalues()).real();
                double min_eigenvalue_topleft = eigenvalues_topleft.minCoeff();
                
                // check which one is lower and keep it
                if (min_eigenvalue_topleft < min_eigenvalue) min_eigenvalue = min_eigenvalue_topleft;
            }


        return min_eigenvalue;

    }





    double EssentialVerification::computePenalisedMatrixMAndMinEigenvalueGeneral(const Vector6 & Lagrange_multipliers, 
                                                                                 const Matrix12 & Q, size_t idx_relaxation)
    {
        /* This function checks the optimality of the solution via the general algorithm 
        in the paper.
        The relaxation is indexed by the PARAM *idx_relaxation* from 1 to 5 (inclusive)
        */

        Matrix12 M;
        M.setZero();

        Matrix12 A, As;
        A.setZero();
        As.setZero();
        unsigned int e1=0; unsigned int e4=1;  unsigned int e7=2;
        unsigned int e2=3; unsigned int e5=4;  unsigned int e8=5;
        unsigned int e3=6; unsigned int e6=7;  unsigned int e9=8;
        unsigned int t1=9; unsigned int t2=10; unsigned int t3=11;

        M = Q;

        Vector7 Lagrange_multipliers_ext;
        // fill the vector
        Lagrange_multipliers_ext.setZero();
        // trivial cases
        if (idx_relaxation == 5)
          Lagrange_multipliers_ext.topRows(6) = Lagrange_multipliers;
        else
        {
          Lagrange_multipliers_ext.topRows(idx_relaxation + 1) = Lagrange_multipliers.topRows(idx_relaxation + 1);
          Lagrange_multipliers_ext.bottomRows(5 - idx_relaxation) = Lagrange_multipliers.bottomRows(5 - idx_relaxation);
        }

        // Create the matrices
        As.setZero();
        A.setZero();
        // t1^2 + t2^2 + t3^2 = 1
        As(t1,t1)=1; As(t2,t2)=1;As(t3,t3)=1;    A = 0.5 * (As + As.transpose());
        M -=  Lagrange_multipliers_ext(0) * A;


        As.setZero();
        A.setZero();
        // e1^2 + e2^2 + e3^2 = t2 ^2 + t3^2
        As(e1,e1)=1; As(e2,e2)=1;As(e3,e3)=1;As(t2, t2)=-1;As(t3, t3)=-1;  A = 0.5 * (As + As.transpose());
        M -=  Lagrange_multipliers_ext(1) * A;


        As.setZero();
        A.setZero();
        // e1*e4 * e2*e5 + e3*e6 = -t1*t2
        As(e1,e4)=1; As(e2,e5)=1;As(e3,e6)=1;As(t1, t2)=1;  A = 0.5 * (As + As.transpose());
        M -=  Lagrange_multipliers_ext(2) * A;

        As.setZero();
        A.setZero();
        // e1*e7 * e2*e8 + e3*e9 = -t1*t3
        As(e1,e7)=1; As(e2,e8)=1;As(e3,e9)=1;As(t1, t3)=1;  A = 0.5 * (As + As.transpose());
        M -=  Lagrange_multipliers_ext(3) * A;

        As.setZero();
        A.setZero();
        // e4^2 + e5^2 + e6^2 = t1 ^2 + t3^2
        As(e4,e4)=1; As(e5,e5)=1;As(e6,e6)=1;As(t1, t1)=-1; As(t3, t3)=-1; A = 0.5 * (As + As.transpose());
        M -=  Lagrange_multipliers_ext(4) * A;

        As.setZero();
        A.setZero();
        // e4*e7 * e5*e8 + e6*e9 = -t2*t3
        As(e4,e7)=1; As(e5,e8)=1;As(e6,e9)=1;As(t2, t3)=1; A = 0.5 * (As + As.transpose());
        M -=  Lagrange_multipliers_ext(5) * A;

        As.setZero();
        A.setZero();
        // e7^2 + e8^2 + e9^2 = t1 ^2 + t2^2
        As(e7,e7)=1; As(e8,e8)=1;As(e9,e9)=1;As(t1, t1)=-1; As(t2, t2)=-1;  A = 0.5 * (As + As.transpose());
        M -=  Lagrange_multipliers_ext(6) * A;


        // Compute min eigenvalue of M

        double min_eigenvalue = -1000;
        Vector12 min_eigenvector;

            Vector3 eigenvalues_bottomrightM = ((M.block<3,3>(9, 9)).eigenvalues()).real();
            min_eigenvalue = eigenvalues_bottomrightM.minCoeff(); 
            // std::cout << "Eigenvalues Mt:\n" << eigenvalues_bottomrightM << std::endl; 
            // std::cout << "Minimum of this vector:\n" << min_eigenvalue << std::endl; 
            // In many cases, this part fails
            if (min_eigenvalue > tau_)
            {
                // Try with the other one

                Vector9 eigenvalues_topleft = ((M.block<9,9>(0, 0)).eigenvalues()).real();
                double min_eigenvalue_topleft = eigenvalues_topleft.minCoeff();
                // check which one is lower and keep it
                // std::cout << "Eigenvalues ME:\n" << eigenvalues_topleft << std::endl; 
                // std::cout << "Minimum of this vector:\n" << min_eigenvalue_topleft << std::endl; 
                if (min_eigenvalue_topleft < min_eigenvalue) min_eigenvalue = min_eigenvalue_topleft;
            }

            
        return min_eigenvalue;

    }









    bool EssentialVerification::checkOptimalitySolution(const Vector6 & Lagrange_multipliers,
                                                        const Matrix12 & Q, double& mu_min, 
                                                        double& dual_gap, double& d_hat,
                                                        size_t idx_relaxation)
    {
         bool is_opt = false;
         // Verify PSD of M
         // check first relaxation
         if (idx_relaxation == 0)
                mu_min = computePenalisedMatrixMAndMinEigenvalue(Lagrange_multipliers, Q);
         else
                mu_min = computePenalisedMatrixMAndMinEigenvalueGeneral(Lagrange_multipliers, Q, idx_relaxation);

         // Compute dual objective
         d_hat = Lagrange_multipliers(0);

         // compute dual gap
         dual_gap = f_hat_ - d_hat;
         // Check conditions for optimality
         // std::cout << "Minimum eigenvalue: " << mu_min << std::endl; 
         if (mu_min > tau_)
            is_opt = true;
        return is_opt;
    }




}  // end of essential namespace
