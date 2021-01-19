cmake_minimum_required(VERSION 3.10)


set(PYBIND11_PYTHON_VERSION ${RELATIVEPOSEPYTHON_VERSION})

set(ROBUST_ESTIMATION ON)

project(relative_pose_python_bindings)

pybind11_add_module(relativeposepython relativeposepython/relativeposepython.cc)

if (ROBUST_ESTIMATION)
        pybind11_add_module(robustrelativeposepython robustrelativeposepython/robustrelativeposepython.cc)
endif()

message(STATUS "Python Interpreter Version: ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}")


if (0)
if (NOT (PYTHON_VERSION_MAJOR EQUAL 2 AND PYTHON_VERSION_MINOR EQUAL 7))
    # Hack: VTK used in PCL might add /usr/include/python2.7 to all targets'
    # INCLUDE_DIRECTORIES properties. We need to remove it.
    get_target_property(RELATIVEPOSEPY_NEW_INCLUDE_DIRS relativeposepython INCLUDE_DIRECTORIES)
    list(FILTER RELATIVEPOSEPY_NEW_INCLUDE_DIRS EXCLUDE REGEX ".*python2.7$")
    set_target_properties(relativeposepython
            PROPERTIES INCLUDE_DIRECTORIES "${RELATIVEPOSEPY_NEW_INCLUDE_DIRS}")
endif ()
endif()

#set_target_properties(relativeposepython
#            PROPERTIES INCLUDE_DIRECTORIES "${RELATIVEPOSEPY_NEW_INCLUDE_DIRS}")
target_link_libraries(relativeposepython PUBLIC Essential)

if (ROBUST_ESTIMATION)
     target_link_libraries(robustrelativeposepython PUBLIC Essential)
endif()

# fix for clang
# see: https://github.com/pybind/pybind11/issues/1818
# if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
#     target_compile_options(relativeposepython PUBLIC -fsized-deallocation)
# endif ()
add_compile_options(-I/usr/include/python2.7 -lpython2.7 )


# make sure to output the build file to relativeposepython folder
SET_TARGET_PROPERTIES(relativeposepython
        PROPERTIES
        PREFIX ""
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/relativeposepython"
        )

if (ROBUST_ESTIMATION)
SET_TARGET_PROPERTIES(robustrelativeposepython
        PROPERTIES
        PREFIX ""
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/robustrelativeposepython"
        )
endif()

# copy package __init__.py file
configure_file(relativeposepython/__init__.py
        ${CMAKE_CURRENT_BINARY_DIR}/relativeposepython/__init__.py
        )

if (ROBUST_ESTIMATION)
configure_file(robustrelativeposepython/__init__.py
        ${CMAKE_CURRENT_BINARY_DIR}/robustrelativeposepython/__init__.py
        )
endif()


# copy setup.py file
configure_file(relativeposepython/setup.py.in
        ${CMAKE_CURRENT_BINARY_DIR}/relativeposepython/setup.py)
        
if (ROBUST_ESTIMATION)       
configure_file(robustrelativeposepython/setup.py.in
        ${CMAKE_CURRENT_BINARY_DIR}/robustrelativeposepython/setup.py)
endif()


file(COPY .
        DESTINATION .
        FILES_MATCHING
        PATTERN *.py)
