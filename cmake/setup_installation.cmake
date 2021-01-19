include(CMakePackageConfigHelpers)
# Include module with fuction 'write_basic_package_version_file'

# Configuration

# Set up install directories. INCLUDE_INSTALL_DIR, LIB_INSTALL_DIR and
# CMAKECONFIG_INSTALL_DIR must not be absolute paths.
set(include_install_dir "include/Essential")
set(LIB_INSTALL_DIR "${LIBRARY_TARGET_NAME}")
set(CMAKECONFIG_INSTALL_DIR "share/${LIBRARY_TARGET_NAME}")
set(RELATIVE_CMAKECONFIG_INSTALL_DIR "share/${LIBRARY_TARGET_NAME}")
set(config_install_dir "lib/cmake/${LIBRARY_TARGET_NAME}-${PROJECT_VERSION}")
set(exported_targets_name "${LIBRARY_TARGET_NAME}Targets")
set(exported_targets_filename "${exported_targets_name}.cmake")
set(export_dirpath "lib/cmake/${LIBRARY_TARGET_NAME}")
set(config_basename "${LIBRARY_TARGET_NAME}Config")
set(config_filename "${config_basename}.cmake")
set(version_filename "${config_basename}Version.cmake")

####################################################
##                 THEIA  INSPIRED                ##
####################################################

# Install the .h files
file(GLOB ${LIBRARY_TARGET_NAME}_HDRS_GLOB ${CMAKE_SOURCE_DIR}/include/*.h) 
install(FILES ${${LIBRARY_TARGET_NAME}_HDRS_GLOB} DESTINATION ${include_install_dir})

install(DIRECTORY /include DESTINATION include/${LIBRARY_TARGET_NAME} FILES_MATCHING PATTERN "*.h")


# This "exports" all targets which have been put into the export set
# "TheiaExport". This means that CMake generates a file with the given
# filename, which can later on be loaded by projects using this package.
# This file contains ADD_LIBRARY(bar IMPORTED) statements for each target
# in the export set, so when loaded later on CMake will create "imported"
# library targets from these, which can be used in many ways in the same way
# as a normal library target created via a normal ADD_LIBRARY().
install(EXPORT ${LIBRARY_TARGET_NAME_EXPORT}
        DESTINATION ${RELATIVE_CMAKECONFIG_INSTALL_DIR} FILE ${exported_targets_filename})

# Figure out the relative path from the installed Config.cmake file to the
# install prefix (which may be at runtime different from the chosen
# CMAKE_INSTALL_PREFIX if under Windows the package was installed anywhere)
# This relative path will be configured into the TheiaConfig.cmake.
file(RELATIVE_PATH INSTALL_ROOT_REL_CONFIG_INSTALL_DIR
     ${CMAKE_INSTALL_PREFIX}/${CMAKECONFIG_INSTALL_DIR} ${CMAKE_INSTALL_PREFIX})

# Create a TheiaConfig.cmake file. <name>Config.cmake files are searched by
# FIND_PACKAGE() automatically. We configure that file so that we can put any
# information we want in it, e.g. version numbers, include directories, etc.
configure_package_config_file(
                                "${CMAKE_SOURCE_DIR}/cmake/${config_filename}.in"
                                "${CMAKE_CURRENT_BINARY_DIR}/${config_filename}"
                                INSTALL_DESTINATION "${config_install_dir}")

# Additionally, when CMake has found a TheiaConfig.cmake, it can check for a
# TheiaConfigVersion.cmake in the same directory when figuring out the version
# of the package when a version has been specified in the FIND_PACKAGE() call,
# e.g. FIND_PACKAGE(Theia [0.5.2] REQUIRED). The version argument is optional.
configure_file("${CMAKE_SOURCE_DIR}/cmake/${version_filename}.in"
               "${CMAKE_CURRENT_BINARY_DIR}/${version_filename}" @ONLY)

# Install these files into the same directory as the generated exports-file,
# we include the FindPackage scripts for libraries whose headers are included
# in the public API of Theia and should thus be present in THEIA_INCLUDE_DIRS.
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${config_filename}"
              "${CMAKE_CURRENT_BINARY_DIR}/${version_filename}"
              "${CMAKE_SOURCE_DIR}/cmake/FindEigen.cmake"
              "${CMAKE_SOURCE_DIR}/cmake/FindOptimization.cmake"
              "${CMAKE_SOURCE_DIR}/cmake/Findgncso.cmake"
              "${CMAKE_SOURCE_DIR}/cmake/FindOpengv.cmake"
              DESTINATION ${CMAKECONFIG_INSTALL_DIR})
  
  
 # uninstall target
if(NOT TARGET uninstall)
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/CMakeUninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake/CMakeUninstall.cmake"
    IMMEDIATE @ONLY)

  add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake/CMakeUninstall.cmake)
endif()

