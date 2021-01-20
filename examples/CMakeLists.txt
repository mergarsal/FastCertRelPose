project(OptCertifier-Examples CXX)


# Find Eigen library
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../cmake/")
find_package(Eigen3 3.3.3 REQUIRED)
if(EIGEN3_FOUND)
message(STATUS "Found Eigen3 library (version ${EIGEN3_VERSION})")
message(STATUS "Eigen3 include directory: ${EIGEN3_INCLUDE_DIR}\n")
else()
message(STATUS "Eigen library not found!")
endif()

# Add the Eigen include directories
include_directories(${EIGEN3_INCLUDE_DIR})

# Fast Essential Matrix Estimation
add_executable(example_essential_matrix ${CMAKE_CURRENT_SOURCE_DIR}/example_essential_matrix.cpp)
target_link_libraries(example_essential_matrix 
                      Essential
                      Optimization)


# GNC-RTNT  Essential Matrix Estimation
add_executable(example_GNC_TNT_essential_matrix ${CMAKE_CURRENT_SOURCE_DIR}/example_GNC_TNT_essential_matrix.cpp) 

target_link_libraries(example_GNC_TNT_essential_matrix 
                Essential 
                gncso
                Optimization)
 
  
