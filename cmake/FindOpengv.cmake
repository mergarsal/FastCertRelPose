find_package(opengv)

# Set standard CMake FindPackage variables if found.
if (opengv_FOUND)
  set(opengv_INCLUDE_DIRS ${opengv_INCLUDE_DIR})
endif (opengv_FOUND)


# Only mark internal variables as advanced if we found opengv, otherwise
# leave it visible in the standard GUI for the user to set manually.
if (opengv_FOUND)
  mark_as_advanced(FORCE opengv_INCLUDE_DIR
    opengv_DIR)
endif (opengv_FOUND)
