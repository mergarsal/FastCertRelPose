# mgarsal 03/02/2020
# PROJECT CONFIGURATION
cmake_minimum_required(VERSION 3.1)
project(opt_certifier_relative_pose LANGUAGES C CXX VERSION 1.0.0)

# set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel.")
set(CMAKE_BUILD_TYPE Release)


set(LIBRARY_TARGET_NAME "Essential")
set(LIBRARY_TARGET_NAME_EXPORT "${LIBRARY_TARGET_NAME}Export")

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON) # We require C++ 14


set(DO_GNC ON)
# build the examples
set(BUILD_${LIBRARY_TARGET_NAME}_EXAMPLE ON)
# build the matlab binding
set(BUILD_${LIBRARY_TARGET_NAME}_MATLAB OFF)
# build the python binding
set(BUILD_${LIBRARY_TARGET_NAME}_PYTHON OFF)


### CMake Cache (build configuration) variables -- these are set interactively in the CMake GUI, and cached in CMakeCache ###
# Directory for built libraries
set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/lib CACHE PATH "The directory in which to place libraries built by this project")
# Directory for built executables
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/bin CACHE PATH "The directory in which to place executables built by this project")

# BUILD CONFIGURATIONS
option(CMAKE_VERBOSE_MAKEFILE "Generate verbose makefiles?" OFF)

set(CODE_PROFILING OFF CACHE BOOL "Turn on code profiling?")
if(${CODE_PROFILING})
        message(STATUS "Turning on code profiling for Essential Matrix Estimation")
endif()

# Add the .cmake files that ship with Eigen3 to the CMake module path (useful for finding other stuff)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/eigen/cmake" CACHE STRING "The CMake module path used for this project")



# FIND EIGEN3
set( ENV{EIGEN3_ROOT_DIR} ${CMAKE_SOURCE_DIR}/eigen)

find_package(Eigen3 3.3 REQUIRED)

if(EIGEN3_FOUND)
        message(STATUS "Found Eigen3 library (version ${EIGEN3_VERSION_STRING})")
        message(STATUS "Eigen3 include directory: ${EIGEN3_INCLUDE_DIR}\n")
else()
        message(STATUS "Eigen library not found!")
endif()


# FIND ADDITIONAL LIBRARIES
# These next operations make use of the .cmake files shipped with Eigen3
find_package(BLAS REQUIRED)

# IMPORT OpenGV library
add_definitions(-march=native)

find_package(opengv REQUIRED)

find_path(opengv_INCLUDE_DIR
  NAMES "opengv/relative_pose/MACentralRelative.hpp"
  PATHS # ${CMAKE_SOURCE_DIR}/opengv/include/opengv
        /usr/include /usr/local/include
        NO_DEFAULT_PATH)
      
if(opengv_INCLUDE_DIR)
  message("-- Found OpenGV headers: ${opengv_INCLUDE_DIR}")
else(opengv_INCLUDE_DIR)
  message(FATAL_ERROR "Can't find OpenGV headers. Try passing -Dopengv_INCLUDE_DIR=...")
endif(opengv_INCLUDE_DIR)


include_directories(${opengv_INCLUDE_DIR})

# Installation of Headers
install(DIRECTORY ${PROJECT_SOURCE_DIR}/opengv/include
  DESTINATION ${CMAKE_INSTALL_PREFIX}
  FILES_MATCHING PATTERN "*.h"
  PATTERN "CMakeFiles" EXCLUDE)
  
  
  
  

add_subdirectory(dependences/Rosen_optimization)
# Installation of Headers
install(  DIRECTORY ${PROJECT_SOURCE_DIR}/dependences/Rosen_optimization/include
          DESTINATION ${CMAKE_INSTALL_PREFIX}
          FILES_MATCHING PATTERN "*.h"
          PATTERN "CMakeFiles" EXCLUDE
)

include_directories(${Optimization_INCLUDE_DIR})


if(DO_GNC)
        find_package(gncso REQUIRED)
endif()




set(${LIBRARY_TARGET_NAME}_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(${LIBRARY_TARGET_NAME}_SOURCE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(${LIBRARY_TARGET_NAME}_EXAMPLES_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/examples)


# Expose the include directories for this project
set(${LIBRARY_TARGET_NAME}_ADD_INCLUDES ${Optimization_INCLUDE_DIR} 
                                        ${EIGEN3_INCLUDE_DIR} 
                                        ${opengv_INCLUDE_DIR} 
                                        ${gncso_INCLUDE_DIR}
                                        ${Optimization_INCLUDE_DIR})
                                
set(${LIBRARY_TARGET_NAME}_CERT_INCLUDES 
        ${${LIBRARY_TARGET_NAME}_INCLUDE_DIRS} 
        ${${LIBRARY_TARGET_NAME}_ADD_INCLUDES})




# Get the set of Essential header and source files
set(${LIBRARY_TARGET_NAME}_HDRS
${${LIBRARY_TARGET_NAME}_INCLUDE_DIRS}/EssentialTypes.h
${${LIBRARY_TARGET_NAME}_INCLUDE_DIRS}/EssentialUtils.h
${${LIBRARY_TARGET_NAME}_INCLUDE_DIRS}/EssentialManifold.h
${${LIBRARY_TARGET_NAME}_INCLUDE_DIRS}/EssentialProblem.h
${${LIBRARY_TARGET_NAME}_INCLUDE_DIRS}/Essential.h
${${LIBRARY_TARGET_NAME}_INCLUDE_DIRS}/CorrespondencesVectors.h
${${LIBRARY_TARGET_NAME}_INCLUDE_DIRS}/generateCorrespondences.h
${${LIBRARY_TARGET_NAME}_INCLUDE_DIRS}/EssentialVerification.h
${${LIBRARY_TARGET_NAME}_INCLUDE_DIRS}/GNCEssential.h
)


set(${LIBRARY_TARGET_NAME}_SRCS
${${LIBRARY_TARGET_NAME}_SOURCE_DIRS}/EssentialUtils.cpp
${${LIBRARY_TARGET_NAME}_SOURCE_DIRS}/EssentialManifold.cpp
${${LIBRARY_TARGET_NAME}_SOURCE_DIRS}/EssentialProblem.cpp
${${LIBRARY_TARGET_NAME}_SOURCE_DIRS}/Essential.cpp
${${LIBRARY_TARGET_NAME}_SOURCE_DIRS}/generateCorrespondences.cpp
${${LIBRARY_TARGET_NAME}_SOURCE_DIRS}/EssentialVerification.cpp
${${LIBRARY_TARGET_NAME}_SOURCE_DIRS}/GNCEssential.cpp
)



# Build the Essential library
add_library(${LIBRARY_TARGET_NAME} 
        ${${LIBRARY_TARGET_NAME}_HDRS} 
        ${${LIBRARY_TARGET_NAME}_SRCS}  )

  
# Set two minimum target properties for the library.
# See https://cmake.org/cmake/help/latest/command/set_target_properties.html

target_include_directories(${LIBRARY_TARGET_NAME} PUBLIC
        # only when building from the source tree
        $<BUILD_INTERFACE:${${LIBRARY_TARGET_NAME}_INCLUDE_DIRS}>
        # only when using the lib from the install path
        $<INSTALL_INTERFACE:include>
        ${${LIBRARY_TARGET_NAME}_ADD_INCLUDES}
)


target_link_libraries(${LIBRARY_TARGET_NAME} 
                        PUBLIC 
                        opengv 
                        ${BLAS_LIBRARIES} 
                        ${M} 
                        ${LAPACK} 
)

if(DO_GNC)
  target_link_libraries(${LIBRARY_TARGET_NAME}
      PUBLIC
        gncso
  )
endif()

if(${CODE_PROFILING})
  set_target_properties(${LIBRARY_TARGET_NAME} PROPERTIES COMPILE_FLAGS "-pg -g" LINK_FLAGS "-pg -g")
endif()

        
set_target_properties(${LIBRARY_TARGET_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)
        
        
install(TARGETS ${LIBRARY_TARGET_NAME}
  EXPORT  ${LIBRARY_TARGET_NAME_EXPORT}
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib${LIB_SUFFIX}
  ARCHIVE DESTINATION lib${LIB_SUFFIX}
  INCLUDES DESTINATION "include"
  PUBLIC_HEADER DESTINATION "include/${LIBRARY_TARGET_NAME}"
 )
  
  

# Build the example executable
IF(BUILD_${LIBRARY_TARGET_NAME}_EXAMPLE)
        message(STATUS "Adding examples to build")
        add_subdirectory(examples)
endif()


# building matlab binding
if(BUILD_${LIBRARY_TARGET_NAME}_MATLAB)
        message(STATUS "Adding matlab binding to build")
        add_subdirectory(matlab)
endif()


# building python binding
if(BUILD_${LIBRARY_TARGET_NAME}_PYTHON)
        # download the pybind11 repo
        configure_file(cmake/pybind11.CMakeLists.txt.in pybind11-download/CMakeLists.txt)
        execute_process(COMMAND "${CMAKE_COMMAND}" -G "${CMAKE_GENERATOR}" .
                    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/pybind11-download")
        execute_process(COMMAND "${CMAKE_COMMAND}" --build .
                    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/pybind11-download")
        add_subdirectory("${CMAKE_BINARY_DIR}/pybind11-src"
                    "${CMAKE_BINARY_DIR}/pybind11-build")
        message(STATUS "Adding python binding to build")
        add_subdirectory(python)
endif()



# Install
include(cmake/setup_installation.cmake)

