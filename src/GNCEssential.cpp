#include "GNCEssential.h"

#include "gncso/Refinement/Base/Concepts.h"
#include "gncso/Refinement/Expansion/ExpansionInliers.h"


using namespace std::chrono;
using namespace Optimization;
using namespace Optimization::Riemannian;
using namespace gncso;


namespace Essential{

    GNCEssentialEstimationResult GNCEssentialClass::getResultGNC(void)
    {

               Matrix9 C;
               if (options_.estimation_verbose)
                                        std::cout << "-----------------\n GNC Essential matrix estimation \n-----------------";

               GNCEssentialEstimationResult result;


               Matrix3 E_8pts;
               Matrix3 matrix_precon;


                //////////////////////////////////////////////////
                if (options_.estimation_verbose)
                    std::cout << " Constructing data matrix C.\n";
                    
                // Get starting timepoint
                auto start_time_init = high_resolution_clock::now();
                // create data matrix
                C = constructDataMatrix(points_);

                auto duration_construction_C = duration_cast<microseconds>(high_resolution_clock::now() - start_time_init);

                auto start_time_bounds = high_resolution_clock::now();
                
                
                //////////////////////////////////////////////////


                // Compute and cache preconditioning matrices
                // Just compute this if required
                if ((options_.use_preconditioning != Preconditioner::Any) || (options_.chosen_initialisation == InitialisationMethod::PTS8_INIT))
                        E_8pts = initialize8pts(C, matrix_precon, options_.use_preconditioning);
              
                auto duration_8pts = duration_cast<microseconds>(high_resolution_clock::now() - start_time_bounds);
                //////////////////////////////////////////////////


                if (options_.estimation_verbose)

                    std::cout << "Initial estimate from linear estimation (8pts): \n" << E_8pts << std::endl;



               if (options_.estimation_verbose)
               {
                     std::cout << "Choosing initial guess: \n" << std::endl;
               }
               // Update initial guess
               Matrix3 E_initial;

               Matrix3 R_s;
               Vector3 t_s;

               auto start_time_init_essential = high_resolution_clock::now();

               switch (options_.chosen_initialisation)
               {
                case InitialisationMethod::RANDOM_INIT:
                    if (options_.estimation_verbose) std::cout << "Random guess... \n";
                    E_initial = projectToEssentialManifold(Matrix3::Random());
                    computeRtfromE(points_, E_initial, R_s, t_s);
                    break;

                case InitialisationMethod::USER_SUPPLIED:
                    if (options_.estimation_verbose) std::cout << "User supplied guess ... \n";
                    E_initial = projectToEssentialManifold(E_initial_);
                    computeRtfromE(points_, E_initial_, R_s, t_s);
                    break;
                    
                case InitialisationMethod::PTS5_INIT:
                    if (options_.estimation_verbose) std::cout << "5 points algorithm init ... \n";
                    E_initial = initialize5pts(points_);
                    computeRtfromE(points_, E_initial, R_s, t_s);
                    break;

                 case InitialisationMethod::PTS7_INIT:
                    if (options_.estimation_verbose) std::cout << "7 points algorithm init ... \n";
                    E_initial = initialize7pts(points_);
                    computeRtfromE(points_, E_initial, R_s, t_s);
                    break;

                case InitialisationMethod::PTS8_INIT:
                    if (options_.estimation_verbose)  std::cout << "8 points guess ... \n";
                    E_initial = E_8pts;
                    computeRtfromE(points_, E_initial, R_s, t_s);
                    break;
                    
                //case InitialisationMethod::EYE_INIT:
                default:
                    if (options_.estimation_verbose) std::cout << "Identity guess... \n";
                    E_initial = projectToEssentialManifold(Matrix3::Identity());
                    computeRtfromE(points_, E_initial, R_s, t_s);
                    break;
               }

              auto duration_init_diff_methods = duration_cast<microseconds>(high_resolution_clock::now() - start_time_init_essential);


               if (options_.estimation_verbose)
               {
                std::cout << "Initial guess rotation:\n" << R_s << std::endl;
                std::cout << "Initial guess translation Tgt:\n" << t_s << std::endl;
               }


               if (options_.estimation_verbose)
                    std::cout << "Recovering R & t from E: \nR : \n" << R_s << "\nt:\n" << t_s << std::endl;


               Matrix34 Rt;
               Rt.block<3, 3>(0, 0) = R_s;
               Rt.block<3, 1>(0, 3) = t_s;
               // initial weights (1)
               weights_t initial_weights(points_.size());
               initial_weights.setOnes();

               // Show initial objective
               double f_initial = 0.5 * vec(E_initial).transpose() * C * vec(E_initial);

               if (options_.estimation_verbose)
                   std::cout << "Initial objective value f_{init} = " << f_initial << std::endl;



	       // Define the problem
	       if (options_.estimation_verbose)         std::cout << "Creating problem " << std::endl;

               /// Define the problem
               EssentialProblem problem(C, options_.use_preconditioning);
               problem.setPointCorrespondences(points_);
               // set precon matrix
               if (options_.use_preconditioning != Preconditioner::None)
                 problem.setMatrixPrecon(matrix_precon);

                // Cache these three matrices for gradient & hessian
                ProblemCachedMatrices problem_matrices;
            
            
            
            
               if (options_.estimation_verbose)         std::cout << "Setting precon\n";
                /// Function handles required by the TNT optimization algorithm
                // Preconditioning operator (optional)
                  std::experimental::optional<Optimization::Riemannian::LinearOperator<Matrix34, Matrix34, ProblemCachedMatrices>> precon;

                  if (options_.use_preconditioning == Preconditioner::None)
                    {
                        precon = std::experimental::nullopt;
                    }
                  else {
                    Optimization::Riemannian::LinearOperator<Matrix34, Matrix34, ProblemCachedMatrices> precon_op =
                        [&problem](const Matrix34 &Y, const Matrix34 &Ydot,
                                   const ProblemCachedMatrices & problem_matrices) {
                          return problem.precondition(Y, Ydot);
                        };
                    precon = precon_op;
                  }





                  if (options_.estimation_verbose)          std::cout << "Setting objective\n";
                  
              // Objective
              Optimization::Objective<Matrix34, double, ProblemCachedMatrices> F =
                  [&problem](const Matrix34 &Y, ProblemCachedMatrices & problem_matrices){
                    return problem.evaluate_objective(Y, problem_matrices);
                  };
                  
                  



                  if (options_.estimation_verbose)         std::cout << "Setting gradient\n";
                  
             /// Gradient
              Optimization::Riemannian::VectorField<Matrix34, Matrix34, ProblemCachedMatrices> grad_F = [&problem](const Matrix34 &Y,
                                                ProblemCachedMatrices & problem_matrices) {
                // Compute and cache Euclidean gradient at the current iterate
                problem_matrices.NablaF_Y = problem.Euclidean_gradient(Y, problem_matrices);
                // Compute Riemannian gradient from Euclidean one
                return problem.Riemannian_gradient(Y, problem_matrices.NablaF_Y);
              };




              if (options_.estimation_verbose)         std::cout << "Setting quadratic model\n";
              
              
              // Local quadratic model constructor
              Optimization::Riemannian::QuadraticModel<Matrix34, Matrix34, ProblemCachedMatrices> QM =
                  [&problem](
                      const Matrix34 &Y, Matrix34 &grad,
                      Optimization::Riemannian::LinearOperator<Matrix34, Matrix34, ProblemCachedMatrices> &HessOp,
                      ProblemCachedMatrices & problem_matrices) {
                    // Compute and cache Euclidean gradient at the current iterate

                    problem_matrices.NablaF_Y  = problem.Euclidean_gradient(Y, problem_matrices);
                    // Compute Riemannian gradient from Euclidean gradient
                    grad = problem.Riemannian_gradient(Y, problem_matrices.NablaF_Y );


                    // Define linear operator for computing Riemannian Hessian-vector
                    // products
                    HessOp = [&problem](const Matrix34 &Y, const Matrix34 &Ydot,
                                        const ProblemCachedMatrices & problem_matrices) {
                                        Matrix34 Hss = problem.Riemannian_Hessian_vector_product(Y, problem_matrices, Ydot);
                                        return Hss;
                    };
                  };



                  // Riemannian metric
                  if (options_.estimation_verbose)         std::cout << "Setting metric\n";
                  // We consider a realization of the product of Stiefel manifolds as an
                  // embedded submanifold of R^{r x dn}; consequently, the induced Riemannian
                  // metric is simply the usual Euclidean inner product
                  Optimization::Riemannian::RiemannianMetric<Matrix34, Matrix34, double, ProblemCachedMatrices>
                      metric = [&problem](const Matrix34 &Y, const Matrix34 &V1, const Matrix34 &V2,
                                          const ProblemCachedMatrices & problem_matrices) {

                        Matrix3 R1, R2;
                        Vector3 t1, t2;
                        R1 = V1.block<3,3>(0,0);
                        R2 = V2.block<3,3>(0,0);

                        t1 = V1.block<3,1>(0,3);
                        t2 = V2.block<3,1>(0,3);

                        return ((R1 * R2.transpose()).trace() + t1.dot(t2));
                      };

                      if (options_.estimation_verbose)       std::cout << "Setting retraction\n";
        
              // Retraction operator
              Optimization::Riemannian::Retraction<Matrix34, Matrix34, ProblemCachedMatrices> retraction =
                  [&problem](const Matrix34 &Y, const Matrix34 &Ydot, const ProblemCachedMatrices & problem_matrices) {
                    return problem.retract(Y, Ydot);
                  };




            if (options_.estimation_verbose)      std::cout << "Setting function for computing residuals\n";

            // Compute residuals (for GNC)
            baseOpt::ComputeResiduals<Matrix34, weights_t, ProblemCachedMatrices> compute_residuals =
            [&problem](const Matrix34 & X, const ProblemCachedMatrices & problem_matrices)
            {
              return problem.computeResiduals(X);
            };
            
            
            
            

            if (options_.estimation_verbose)      std::cout << "Setting function for updating the weights inside the object\n";

            std::experimental::optional<baseOpt::UpdateParamsInProblem<Matrix34, double, weights_t, ProblemCachedMatrices>> update_params_in_problem_op;
            baseOpt::UpdateParamsInProblem<Matrix34, double, weights_t, ProblemCachedMatrices> update_params_in_problem = [&problem](const Matrix34 & X, const weights_t& new_weights, const ProblemCachedMatrices & problem_matrices)
            {
              return problem.updateWeights(new_weights);
            };

            update_params_in_problem_op = update_params_in_problem;



            if (options_.estimation_verbose)      std::cout << "Setting function for initializing the solver\n";

            std::experimental::optional<baseOpt::InitializeVariableX<Matrix34, double, weights_t, ProblemCachedMatrices>>  initialize_solver_8pt_op;
            baseOpt::InitializeVariableX<Matrix34, double, weights_t, ProblemCachedMatrices>  initialize_solver_8pt = [&problem](const Matrix34 & X, const weights_t& weights, const ProblemCachedMatrices & problem_matrices)
            {
              return problem.initializeSolver8pts();
            };
            initialize_solver_8pt_op = initialize_solver_8pt;



              // Stop timer
              auto stop_init = high_resolution_clock::now();

              auto duration_init = duration_cast<microseconds>( stop_init- start_time_init);




              if (options_.estimation_verbose)     std::cout << "Creating struct with params\n";
              // set up params for solver
              gncso::GNCSmoothParams<double> params;
              params.gradient_tolerance = options_.tol_grad_norm;
              params.relative_decrease_tolerance = options_.rel_func_decrease_tol;
              params.max_iterations = options_.max_RTR_iterations;
              params.max_TPCG_iterations = options_.max_tCG_iterations;
              params.preconditioned_gradient_tolerance = options_.preconditioned_grad_norm_tol;
              params.stepsize_tolerance = options_.stepsize_tol;

              // Set up parameters for GNC
              params.max_outer_iterations = options_.max_outer_iterations;
              params.max_inner_iterations = options_.max_inner_iterations;
              params.cost_diff_threshold = options_.cost_diff_threshold;
              params.mu_threshold = options_.mu_threshold;
              params.max_res_tol_sq = options_.max_res_tol_sq;
              params.inliers_threshold = options_.inliers_threshold;
              params.gnc_factor = options_.gnc_factor;
              params.log_outer_iters = options_.GNClog_iterates;

              params.nr_min_points = options_.nr_min_points;

              // verbose
              params.verbose = options_.estimation_verbose;
              params.GNC_verbose = options_.GNC_verbose;
          

              /** An optional user-supplied function that can be used to instrument/monitor
               * the performance of the internal Riemannian truncated-Newton trust-region
               * optimization algorithm as it runs. */
              // std::experimental::optional<EssentialTNTUserFunction> user_fcn;

              /// Run optimization!
              // GNC
              result.certifier_status = EssentialEstimationStatus::RS_ITER_LIMIT;


              if (options_.estimation_verbose)        std::cout << "Calling GNC function\n";

              GNCResult<Matrix34, weights_t, double> GNCResults;

              auto start_opt = high_resolution_clock::now();
              
              switch (options_.gnc_robust)
              {
               case GNCRobustFunction::TLS:
                   if (options_.estimation_verbose) std::cout << "Calling TLS-GNC... \n";
                   GNCResults = TLSGNCSmooth<Matrix34, weights_t, Matrix34, double, ProblemCachedMatrices>(
                   F, QM, metric, retraction, compute_residuals, Rt, initial_weights, problem_matrices, update_params_in_problem_op, initialize_solver_8pt_op, precon, params);
                   break;
                   
               case GNCRobustFunction::TEMP:
                   if (options_.estimation_verbose) std::cout << "Calling TEMPERATURE ... \n";
                   GNCResults = TemperatureSmooth<Matrix34, weights_t, Matrix34, double, ProblemCachedMatrices>(
                   F, QM, metric, retraction, compute_residuals, Rt, initial_weights, problem_matrices, update_params_in_problem_op, initialize_solver_8pt_op, precon, params);
                   break;
                   
               case GNCRobustFunction::WELSCH:
                   if (options_.estimation_verbose) std::cout << "Calling Welsch ... \n";
                   GNCResults = WelschGNCSmooth<Matrix34, weights_t, Matrix34, double, ProblemCachedMatrices>(
                   F, QM, metric, retraction, compute_residuals, Rt, initial_weights, problem_matrices, update_params_in_problem_op, initialize_solver_8pt_op, precon, params);
                   break;
                   
               case GNCRobustFunction::TUKEY:
                   if (options_.estimation_verbose) std::cout << "Calling Tukey ... \n";
                   GNCResults = TukeyGNCSmooth<Matrix34, weights_t, Matrix34, double, ProblemCachedMatrices>(
                   F, QM, metric, retraction, compute_residuals, Rt, initial_weights, problem_matrices, update_params_in_problem_op, initialize_solver_8pt_op, precon, params);
                   break;
                   
               default: // GNCRobustFunction::GM
                   if (options_.estimation_verbose) std::cout << "Calling GM-GNC... \n";
                   GNCResults = GMGNCSmooth<Matrix34, weights_t, Matrix34, double, ProblemCachedMatrices>(
                   F, QM, metric, retraction, compute_residuals, Rt, initial_weights, problem_matrices, update_params_in_problem_op, initialize_solver_8pt_op, precon, params);
               break;
              }

              auto stop_opt = high_resolution_clock::now();
              auto duration_opt = duration_cast<microseconds>(stop_opt - start_opt);


               // Extract the results
               if (options_.estimation_verbose)
               {
                std::cout << "Final point SO(3) x S(2): " << GNCResults.x << std::endl;
                std::cout << "Final objective: " << GNCResults.f << std::endl;
                std::cout << "Final mu: " << GNCResults.mu << std::endl;
                std::cout << "Final outer iteration: " << GNCResults.nr_outer_iterations << std::endl;
               }

	           Matrix34 Rs_opt = GNCResults.x;
	           Matrix3 R_opt = Rs_opt.block<3,3>(0,0);
	           Vector3 t_opt = Rs_opt.block<3,1>(0,3);


	           Matrix3 E_opt = computeEfromRt(R_opt, t_opt);

	           double f_hat = GNCResults.f;
	           double mu_min, dual_gap, d_hat;
	           Matrix12 Q;


                   if (options_.estimation_verbose)       std::cout << "Cerfication...\n";
	           /// Verification
	           // we need to compute again C here
                   problem.updateWeights(GNCResults.set_inliers);

	           EssentialVerification verify_solution(problem.getDataMatrixC(), E_opt, t_opt, f_hat, options_.eig_min_tol, options_.dual_gap_tol);

                  
                   // output from the verification
                   Vector6 Lagrange_multipliers; 
                   bool is_opt_bool; 
                   size_t idx_relaxation = 0;
                   
                   // total times 
                   double lagrange_total = 0.0; 
                   double verification_total = 0.0; 
                   
                   
                   if (options_.use_all_relaxations == false)
                   {
                   
                           if (options_.use_idx_relaxation == 0)
                           {
	                           // 1. Compute Lagrange multipliers

                                   if (options_.estimation_verbose) std::cout << "Checking only one certifier\n"; 
                                   
                                   auto start_lagrange_reduced = high_resolution_clock::now();
                                   
	                           Lagrange_multipliers = verify_solution.computeLagMult(Q);
	                           
	                           auto duration_lagrange_reduced = duration_cast<microseconds>(high_resolution_clock::now() - start_lagrange_reduced);
	                           
	                           lagrange_total = duration_lagrange_reduced.count(); 
	                          	               


	                           // 2. Check optimality (PSD)
                                   auto start_verification_reduced = high_resolution_clock::now();
                                   
                                   is_opt_bool = verify_solution.checkOptimalitySolution(Lagrange_multipliers, Q, mu_min, dual_gap, d_hat);
                                   
                                   auto duration_verification_reduced = duration_cast<microseconds>(high_resolution_clock::now() - start_verification_reduced);
                                   
                                   
                                   verification_total = duration_verification_reduced.count();
                                   
                                   if (is_opt_bool)     idx_relaxation = 0;   // save relaxation 0 (original)
                                   
                            }
                           
                            else 
                            {
                                if (options_.estimation_verbose)  std::cout << "Checking certifier for relaxation " << 
                                                                                        options_.use_idx_relaxation << "\n";
                                  // for sanity 
                                  if (options_.use_idx_relaxation > 5) 
                                        options_.use_idx_relaxation = 5; 
                        
                                  idx_relaxation = options_.use_idx_relaxation; 
                                 
                                  auto start_lagrange_reduced = high_resolution_clock::now();
                                
                                  Lagrange_multipliers = verify_solution.computeLagMultGeneral(Q, idx_relaxation);
                                  
                                  auto duration_lagrange_reduced = duration_cast<microseconds>(high_resolution_clock::now() - start_lagrange_reduced);
                                  
	                          lagrange_total = duration_lagrange_reduced.count(); 
	                          
	                           
                                  auto start_verification_reduced = high_resolution_clock::now();
                                  
                                  is_opt_bool = verify_solution.checkOptimalitySolution(Lagrange_multipliers, Q, mu_min, dual_gap, d_hat, idx_relaxation);
                                  
                                  auto duration_verification_reduced = duration_cast<microseconds>(high_resolution_clock::now() - start_verification_reduced);
                                  
                                  verification_total = duration_verification_reduced.count();
                            
                            }
                    }
                    // if use all the relaxations
                    else
                    {
                       /** Checking additional certifiers **/
                       
                       
                       if (options_.estimation_verbose)  std::cout << "Checking multiple multipliers\n"; 
                       
                       // if the user asks for all the relaxations
                       if (options_.use_idx_relaxation > 5) 
                        options_.use_idx_relaxation = 5; 
                        
                       for (idx_relaxation = 0; idx_relaxation <= options_.use_idx_relaxation; idx_relaxation++)
                       {
                          auto start_lagrange_reduced = high_resolution_clock::now();
                          Lagrange_multipliers = verify_solution.computeLagMultGeneral(Q, idx_relaxation);
                          auto duration_lagrange_reduced = duration_cast<microseconds>(high_resolution_clock::now() - start_lagrange_reduced);
	                  lagrange_total += duration_lagrange_reduced.count(); 
	                   
                          auto start_verification_reduced = high_resolution_clock::now();
                          is_opt_bool = verify_solution.checkOptimalitySolution(Lagrange_multipliers, Q, mu_min, dual_gap, d_hat, idx_relaxation);
                          auto duration_verification_reduced = duration_cast<microseconds>(high_resolution_clock::now() - start_verification_reduced);
                          verification_total += duration_verification_reduced.count();
                          
                          if (is_opt_bool == 1)
                                break;  // keep the last relaxation
                       }
                       
                       
                     }  // end- select relaxation for the certifier
                     
                    if (is_opt_bool == 0)       idx_relaxation = 10;  // we could not certify the solution 
                    

                 auto duration_total = duration_cast<microseconds>(high_resolution_clock::now() - start_time_init);
             
                if (options_.estimation_verbose)   std::cout << "Tight relaxation with index: " << idx_relaxation << std::endl;
            
 
              
                      // extract results 
                      Rs_opt = GNCResults.x;
	              R_opt = Rs_opt.block<3,3>(0,0);
	              t_opt = Rs_opt.block<3,1>(0,3);
                      E_opt = computeEfromRt(R_opt, t_opt);                      
	           
                      // save GNC
                      result.weights = GNCResults.weights;
                      result.mu = GNCResults.mu;
                      // inliers
                      result.set_inliers = GNCResults.set_inliers;
                      result.elapsed_time = GNCResults.elapsed_time;
                      result.GNCstatus = GNCResults.GNCstatus;
                      result.nr_outer_iterations = GNCResults.nr_outer_iterations;
                      result.valid_estimation = GNCResults.valid_estimation; 
                      
                      // save results for each outer iteration
                      
                      if (options_.GNClog_iterates)
                      {
                        for (size_t id_out = 0; id_out < GNCResults.intermediate_outer_results.size(); id_out++)
                                result.intermediate_outer_results.push_back(GNCResults.intermediate_outer_results[id_out]);
                        
                      }
                      
                      


          	       //  assign all the results to this struct
          	       result.idx_relaxation = idx_relaxation; 
      	               result.f_hat = f_hat;
      	               result.d_hat = d_hat;
      	               result.dual_gap = dual_gap;
      	               result.gradnorm = problem.Riemannian_gradient(Rs_opt).norm();
      	               result.mu_min = mu_min;
      	               result.lagrange_multipliers = Lagrange_multipliers;
                       result.R_opt = R_opt;
                       result.t_opt = t_opt;
                       result.E_opt = E_opt;
                       
                       // Save time
                       result.elapsed_init_time = duration_init.count();
                       result.elapsed_C_time = duration_construction_C.count();  // in microsecs
                       result.elapsed_8pt_time = duration_8pts.count();  // in microsecs
                       result.elapsed_iterative_time = duration_opt.count();  // in microsecs
                       result.elapsed_lagrange_time = lagrange_total;
                       result.elapsed_certifier_time = verification_total;
                       result.elapsed_estimation_time = duration_total.count();
                       result.elapsed_time_methods = duration_init_diff_methods.count();
                       
                       result.certifier_status = EssentialEstimationStatus::RS_ITER_LIMIT;



	                    if (is_opt_bool == true)  result.certifier_status = EssentialEstimationStatus::GLOBAL_OPT;
	                    else result.certifier_status = EssentialEstimationStatus::NO_CERTIFIED;



               return result;

    }; // end of fcn getResult
    
    
    

void GNCEssentialClass::printResultGNC(GNCEssentialEstimationResult& my_result)
{
    // print data

      std::cout << "Data from estimation\n--------------------\n";
      std::cout << "f_hat = " << my_result.f_hat << std::endl;
      std::cout << "d_hat = " << my_result.d_hat << std::endl;

      std::cout << "dual gap: " << my_result.dual_gap << std::endl;
      std::cout << "min. eigenvalue of M = " << my_result.mu_min << std::endl;


      std::cout << "######## Times [in microseconds]:\n";
      std::cout << "Data Matric C construction: " << my_result.elapsed_C_time << std::endl;
      std::cout << "-----\n Total init: " << my_result.elapsed_init_time << std::endl;

      std::cout << "Iterative Method: " << my_result.elapsed_iterative_time << std::endl;

      std::cout << "Computing Lagrange multipliers: " << my_result.elapsed_lagrange_time << std::endl;
      std::cout << "Certifying optimality: " << my_result.elapsed_certifier_time << std::endl;
      std::cout << "Verification: " << my_result.elapsed_certifier_time + my_result.elapsed_lagrange_time << std::endl;
      std::cout << "---------------------\n";
      std::cout << "Total time: " << my_result.elapsed_estimation_time << std::endl << std::endl;

      std::cout << "\n Recovered R:\n" << my_result.R_opt << std::endl;
      std::cout << "\n Recovered t:\n" << my_result.t_opt << std::endl;


      std::cout << "--------\nGNC\n--------\n";
      std::cout << "Final mu: " << my_result.mu << std::endl;

      std::cout << "Total elapsed time: " << my_result.elapsed_time << std::endl;

      std::cout << "Number outer iterations: " << my_result.nr_outer_iterations << std::endl;

      std::cout << "\n Status estimation:\n 0-> GLOBAL OPT \n 1->INCONCLUSIVE \n 2->MAX ITER\n-1-> UNKNOWN STATUS\n\n";

      int is_opt = -1;
      switch(my_result.certifier_status)
      {
        case EssentialEstimationStatus::RS_ITER_LIMIT:
            is_opt = 2;
            break;
        case EssentialEstimationStatus::GLOBAL_OPT:
            is_opt = 0;
            break;
        case EssentialEstimationStatus::NO_CERTIFIED:
            is_opt = 1;
            break;
        default:
            is_opt = -1;
      }

      std::cout << "STATUS: " << is_opt << std::endl;



};  //end of print fcn

}  // end of essential namespace
