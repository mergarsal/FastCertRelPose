find_package(Optimization)

# Set standard CMake FindPackage variables if found.
if (Optimization_FOUND)
  set(Optimization_INCLUDE_DIRS ${Optimization_INCLUDE_DIR})
endif (Optimization_FOUND)


# Only mark internal variables as advanced if we found Eigen, otherwise
# leave it visible in the standard GUI for the user to set manually.
if (Optimization_FOUND)
  mark_as_advanced(FORCE Optimization_INCLUDE_DIR
    Optimization_DIR) 
endif (Optimization_FOUND)
