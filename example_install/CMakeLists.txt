project(Essential-Example-Inst CXX)

cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 14) # We require C++ 14 or later
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_BUILD_TYPE Release)


# Find Eigen library
find_package(Eigen3 3.3.3 REQUIRED)
if(EIGEN3_FOUND)
        message(STATUS "Found Eigen3 library (version ${EIGEN3_VERSION})")
        message(STATUS "Eigen3 include directory: ${EIGEN3_INCLUDE_DIR}\n")
else()
        message(STATUS "Eigen library not found!")
endif()

# Add the Eigen include directories
include_directories(${EIGEN3_INCLUDE_DIR})



  
find_package(Optimization REQUIRED)

find_package(gncso REQUIRED)

find_package(Essential REQUIRED)



# GNC-RTNT  Essential Matrix Estimation
add_executable(example_GNC_TNT_essential_matrix ${CMAKE_CURRENT_SOURCE_DIR}/example_GNC_TNT_essential_matrix.cpp)

# order is important!
target_link_libraries(example_GNC_TNT_essential_matrix 
                Essential gncso )
