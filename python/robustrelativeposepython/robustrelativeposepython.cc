
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include "GNCEssential.h"

namespace py = pybind11;

/**
 * Python interface with pybind11
 */
PYBIND11_MODULE(robustrelativeposepython, m) {
  m.doc() = "Python binding for Robust (Central) Relative Pose";

  

  // Python bound for Essential::EssentialClass
  py::class_<Essential::GNCEssentialClass> essential(m, "GNCEssentialClass");


  // Python bound for Essential::EssentialClass functions
  essential.def(py::init<>())
        .def(py::init< Eigen::Matrix<double, 3, Eigen::Dynamic>&,
                                      Eigen::Matrix<double, 3, Eigen::Dynamic>&, 
                                      Eigen::Matrix<double, 1, Eigen::Dynamic>&, 
                                      const Essential::GNCEssentialEstimationOptions&, 
                                      const Eigen::Matrix<double, 3, 3>&> () )
      .def("getResultGNC", &Essential::GNCEssentialClass::getResultGNC);
      
      

      // Python bound for Essential::InitialisationMethod
      /* Initialization */
      /*  enum class InitialisationMethod {
                            PTS8_INIT = 0,
                            PTS7_INIT,
                            RANDOM_INIT,
                            EYE_INIT,
                            PTS5_INIT,
                            USER_SUPPLIED,
                        };
      */
      py::enum_<Essential::InitialisationMethod>(
      essential, "InitialisationMethod")
      .value("PTS8_INIT", Essential::InitialisationMethod::PTS8_INIT)
      .value("PTS7_INIT", Essential::InitialisationMethod::PTS7_INIT)
      .value("PTS5_INIT", Essential::InitialisationMethod::PTS5_INIT)
      .value("RANDOM_INIT", Essential::InitialisationMethod::RANDOM_INIT)
      .value("EYE_INIT", Essential::InitialisationMethod::EYE_INIT)
      .value("USER_SUPPLIED", Essential::InitialisationMethod::USER_SUPPLIED);
      
      
      /* GNC robust */
      /*
        enum class GNCRobustFunction {
                         None = 0,
                         TLS,
                         GM,
                         TEMP,
                         WELSCH,
                         TUKEY,
                         Any
                     };
      */
      // Python bound for Essential::GNCRobustFunction
      py::enum_<Essential::GNCRobustFunction>(
      essential, "GNCRobustFunction")
      .value("None", Essential::GNCRobustFunction::None)
      .value("TLS", Essential::GNCRobustFunction::TLS)
      .value("GM", Essential::GNCRobustFunction::GM)
      .value("TEMP", Essential::GNCRobustFunction::TEMP)
      .value("WELSCH", Essential::GNCRobustFunction::WELSCH)
      .value("TUKEY", Essential::GNCRobustFunction::TUKEY)
      .value("Any", Essential::GNCRobustFunction::Any);
      

     /* Preconditioner */
     /*
          enum class Preconditioner {
                  None= 0,
                  Max_eigenvalue,
                  Dominant_eigenvalues, 
                  N,
                  Any,
               };
      */
      // Python bound for Essential::GNCRobustFunction
      py::enum_<Essential::Preconditioner>(
      essential, "Preconditioner")
      .value("None", Essential::Preconditioner::None)
      .value("Max_eigenvalue", Essential::Preconditioner::Max_eigenvalue)
      .value("Dominant_eigenvalues", Essential::Preconditioner::Dominant_eigenvalues)
      .value("N", Essential::Preconditioner::N)
      .value("Any", Essential::Preconditioner::Any);
      
      
      
  // Python bound for Essential::EssentialEstimationOptions
  py::class_<Essential::GNCEssentialEstimationOptions>(essential, "Options")
      .def(py::init<>())
      .def_readwrite("estimation_verbose", &Essential::GNCEssentialEstimationOptions::estimation_verbose)
      .def_readwrite("max_res_tol_sq", &Essential::GNCEssentialEstimationOptions::max_res_tol_sq)
      .def_readwrite("gnc_factor", &Essential::GNCEssentialEstimationOptions::gnc_factor)
      .def_readwrite("max_inner_iterations", &Essential::GNCEssentialEstimationOptions::max_inner_iterations)
      .def_readwrite("max_outer_iterations", &Essential::GNCEssentialEstimationOptions::max_outer_iterations)
      .def_readwrite("cost_diff_threshold", &Essential::GNCEssentialEstimationOptions::cost_diff_threshold)
      .def_readwrite("inliers_threshold", &Essential::GNCEssentialEstimationOptions::inliers_threshold)
      .def_readwrite("nr_min_points", &Essential::GNCEssentialEstimationOptions::nr_min_points)
      .def_readwrite("GNC_verbose", &Essential::GNCEssentialEstimationOptions::GNC_verbose)
      .def_readwrite("gnc_robust", &Essential::GNCEssentialEstimationOptions::gnc_robust)
      .def_readwrite("chosen_initialisation", &Essential::GNCEssentialEstimationOptions::chosen_initialisation)  
      .def_readwrite("use_preconditioning", &Essential::GNCEssentialEstimationOptions::use_preconditioning)
      .def_readwrite("use_idx_relaxation", &Essential::EssentialEstimationOptions::use_idx_relaxation);
      
      

      
      
      
      
      
  // Python bound for Essential::EssentialEstimationResult
  py::class_<Essential::GNCEssentialEstimationResult>(essential, "EssentialResult")
      .def(py::init<>())
      .def_readwrite("E_opt", &Essential::GNCEssentialEstimationResult::E_opt)
      .def_readwrite("R_opt", &Essential::GNCEssentialEstimationResult::R_opt)
      .def_readwrite("t_opt", &Essential::GNCEssentialEstimationResult::t_opt)
      .def_readwrite("set_inliers", &Essential::GNCEssentialEstimationResult::set_inliers)
      .def_readwrite("elapsed_estimation_time", &Essential::GNCEssentialEstimationResult::elapsed_estimation_time)
      .def("getStatus", [](const Essential::GNCEssentialEstimationResult& r)
      {
      int is_opt;
      
      
      switch(r.certifier_status)
      {
        case Essential::EssentialEstimationStatus::RS_ITER_LIMIT:
            is_opt = 2;
            break;
        case Essential::EssentialEstimationStatus::GLOBAL_OPT:
            is_opt = 1;
            break;
        case Essential::EssentialEstimationStatus::NO_CERTIFIED:
            is_opt = 0;
            break;
        default:
            is_opt = -1;
      }
      
      
      return is_opt;
      });
      
      
}
