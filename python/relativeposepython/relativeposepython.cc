
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include "Essential.h"

namespace py = pybind11;

/**
 * Python interface with pybind11
 */
PYBIND11_MODULE(relativeposepython, m) {
  m.doc() = "Python binding for (Central) Relative Pose";

  

  // Python bound for Essential::EssentialClass
  py::class_<Essential::EssentialClass> essential(m, "EssentialClass");


  // Python bound for Essential::EssentialClass functions
  essential.def(py::init<>())
        .def(py::init< Eigen::Matrix<double, 3, Eigen::Dynamic>&,
                                      Eigen::Matrix<double, 3, Eigen::Dynamic>&, 
                                      Eigen::Matrix<double, 1, Eigen::Dynamic>&, 
                                      const Essential::EssentialEstimationOptions&, 
                                      const Eigen::Matrix<double, 3, 3>&> () )
      .def("getResults", &Essential::EssentialClass::getResults);
      
      
      
      
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
  py::class_<Essential::EssentialEstimationOptions>(essential, "Options")
      .def(py::init<>())
      .def_readwrite("estimation_verbose", &Essential::EssentialEstimationOptions::estimation_verbose)
      .def_readwrite("chosen_initialisation", &Essential::EssentialEstimationOptions::chosen_initialisation)
      .def_readwrite("use_preconditioning", &Essential::EssentialEstimationOptions::use_preconditioning)
      .def_readwrite("use_idx_relaxation", &Essential::EssentialEstimationOptions::use_idx_relaxation);
      
      
  // Python bound for Essential::EssentialEstimationResult
  py::class_<Essential::EssentialEstimationResult>(essential, "EssentialResult")
      .def(py::init<>())
      .def_readwrite("E_opt", &Essential::EssentialEstimationResult::E_opt)
      .def_readwrite("R_opt", &Essential::EssentialEstimationResult::R_opt)
      .def_readwrite("t_opt", &Essential::EssentialEstimationResult::t_opt)
      //.def_readwrite("certifier_status", &Essential::EssentialEstimationResult::certifier_status)
      .def_readwrite("elapsed_estimation_time", &Essential::EssentialEstimationResult::elapsed_estimation_time)
      .def("getStatus", [](const Essential::EssentialEstimationResult& r)
      {
      int is_opt;
      
      
      switch(r.certifier_status)
      {
        case Essential::EssentialEstimationStatus::RS_ITER_LIMIT:
            is_opt = 2;
            break;
        case Essential::EssentialEstimationStatus::GLOBAL_OPT:
            is_opt = 0;
            break;
        case Essential::EssentialEstimationStatus::NO_CERTIFIED:
            is_opt = 1;
            break;
        default:
            is_opt = -1;
      }
      
      
      return is_opt;
      });
      
      
}
