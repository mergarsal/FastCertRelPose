find_package(gncso)

# Set standard CMake FindPackage variables if found.
if (gncso_FOUND)
  set(gncso_INCLUDE_DIRS ${gncso_INCLUDE_DIR})
  message(STATUS "Found GNCSO in ${GNCSO_INCLUDE_DIRS}")
else()
        message(STATUS "GNCSO not found!")
endif (gncso_FOUND)


# Only mark internal variables as advanced if we found Eigen, otherwise
# leave it visible in the standard GUI for the user to set manually.
if (gncso_FOUND)
  mark_as_advanced(FORCE gncso_INCLUDE_DIR
    gncso_DIR)
endif (gncso_FOUND)
