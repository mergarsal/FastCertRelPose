cmake_minimum_required(VERSION 3.5)

# Following TEASER++
message(STATUS "-- Building matlab binding for essential matrix estimation!")


project(essential_matrix_estimation_matlab_bindings)

find_package(Matlab COMPONENTS MX_LIBRARY)

if (Matlab_FOUND)
    message(STATUS "MATLAB root directory found: ${Matlab_ROOT_DIR}.")

    # no GNC
    matlab_add_mex(
            NAME EssentialMatrixMex
            OUTPUT_NAME EssentialMatrixMex
            SRC EssentialMatrixMex.cc
            LINK_TO Eigen3::Eigen Essential  
    )
    set_target_properties(EssentialMatrixMex PROPERTIES COMPILE_FLAGS "-fvisibility=default")
    
    
    
    # GNC
    matlab_add_mex(
            NAME GNCEssentialMatrixMex
            OUTPUT_NAME GNCEssentialMatrixMex
            SRC GNCEssentialMatrixMex.cc
            LINK_TO Eigen3::Eigen Essential 
    )
    set_target_properties(GNCEssentialMatrixMex PROPERTIES COMPILE_FLAGS "-fvisibility=default")
    
    
    
    
    # copy MATLAB .m files to binary directory
    file(COPY .
            DESTINATION .
            FILES_MATCHING
            PATTERN *.m)
else ()
    message(WARNING "MATLAB root directory not found. Failed to build MATLAB bindings.")
    
endif ()
