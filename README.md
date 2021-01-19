1. Create folder build
``` mkdir build ```
2. Go to build and compile
``` cd build & cmake .. ```
3. 
``` make ```
4.``` cd bin & ./essential_matrix ```


# Known-issues about bindings
## Matlab
Current mex files generator only spport up to gcc4.9. 

## Python:
I dont remember (sth??)





# How to integrate the ``essential matrix estimation ``
# into your project  

## Part 1: Install our library 
0. Create a build folder inside the project folder.

Inside the ``build`` folder, do:
1. cmake .. & make
2. sudo make install

By now you should have everything set up with our library 

## Part 2: Use it in your project
1. Include in your CMakeLists.txt the followig lines: 
``
ADD_DEFINITIONS ( -DPQXX_HIDE_EXP_OPTIONAL )
add_subdirectory(dependences/opt_certifier_relative_pose/C++/Rosen_optimization)
add_subdirectory(dependences/opt_certifier_relative_pose/C++/GNCSO)

# IMPORT OpenGV library
add_definitions(-march=native)
find_package(opengv REQUIRED)

## add essential matrix estimation lib
find_package(Essential REQUIRED)
``

2. [NOT SURE] Add the lines 
```
# Expose the include directories for this project
set(${LIBRARY_TARGET_NAME}_ADD_INCLUDES ${ESSENTIAL_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR} )
```

3. And add to ``target_link_libraries`` the item: ``Essential``

Done!!
