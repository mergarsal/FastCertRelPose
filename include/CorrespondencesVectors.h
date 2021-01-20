#pragma once

#include <iostream>
#include <eigen3/Eigen/Dense>

#include "EssentialTypes.h"

namespace Essential{

        /** A simple struct that contains the elements of a corresponding fatures as bearing (unit) vectors */
        struct CorrespondingFeatures {
         EIGEN_MAKE_ALIGNED_OPERATOR_NEW 
          Vector3 bearing_vector_0, bearing_vector_1;

          double weight_;
         /** Simple default constructor; does nothing */
          CorrespondingFeatures() {}

          CorrespondingFeatures(const Vector3 & bearing_vector_0,
                                const Vector3 & bearing_vector_1,
                                double weight_match = 1.0)
                                : weight_(weight_match), bearing_vector_0(bearing_vector_0),
                                bearing_vector_1(bearing_vector_1) {}

        }; // end of CorrespondingFeatures struct

       // Define types
       typedef std::vector<CorrespondingFeatures> bearing_vectors_t;
       typedef Eigen::VectorXd weights_t;

}  // end of Esssential namespace

