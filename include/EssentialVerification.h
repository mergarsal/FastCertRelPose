#pragma once

#include "EssentialTypes.h"
#include "EssentialUtils.h"
#include <Eigen/Eigenvalues>


namespace Essential{
class EssentialVerification{

    public:
    
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        /* Default constructor*/
        EssentialVerification(void);

        /* Actual constructor */

        EssentialVerification(const Matrix9 & C, const Matrix3 & E,
                        const Vector3 & t, const double f_hat,
                        const double tau = -0.02,
                        const double epsilon_dual = 1e-14):
        C_(C), E_(E), t_(t), f_hat_(f_hat), tau_(tau), epsilon_dual_(epsilon_dual) {};


        /* Default destructor */
        ~EssentialVerification(void){};


        /* Compute lagrange multipliers */
        Vector6 computeLagMult(Matrix12 & Q);
        
        // compute batch of lagrange multipliers (all the relaxations)
        Vector6 computeLagMultGeneral(Matrix12 & Q, size_t idx_relaxation);
        
       /* Certify optimality */
       // NOTE: include all the relaxations
       bool checkOptimalitySolution(const Vector6 & Lagrange_multipliers,
                                    const Matrix12 & Q, double& mu_min, 
                                    double& dual_gap, double& d_hat,
                                    size_t idx_relaxation = 0);


       /* Compute penalised matrix M */
       double computePenalisedMatrixMAndMinEigenvalue(const Vector6 & Lagrange_multipliers, const Matrix12 & Q);
       
       /* compute Hessian of Lagrangian for each relaxation */
       double computePenalisedMatrixMAndMinEigenvalueGeneral(const Vector6 & Lagrange_multipliers, 
                                                        const Matrix12 & Q, size_t idx_relaxation);



    private:
        Vector6 Lagrange_multipliers_;
        // Min eigenvalue
        double mu_min_;
        // Dual gap
        double dual_gap_;

       // Value for the primal & dual objectives
        double f_hat_;
        double dual_hat_;


        // Thresholds
        double tau_;
        double epsilon_dual_;

        // The problem matrices
        Matrix9 C_;
        Matrix3 E_;
        Vector3 t_;




}; // end of Essential verification clas
}  // end of essential namespace
