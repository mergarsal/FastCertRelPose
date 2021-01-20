# Fast and Robust Relative Pose Estimation for Calibrated Cameras

This repository contains the code 
for the relative pose estimation 
between two central and calibrated cameras 
for the [paper](ADD ARXIV)[1]. 

**Authors:** [Mercedes Garcia-Salguero](http://mapir.uma.es/mapirwebsite/index.php/people/290), [Javier Gonzalez-Jimenez](http://mapir.isa.uma.es/mapirwebsite/index.php/people/95-javier-gonzalez-jimenez)


**License:** [GPLv3](https://raw.githubusercontent.com/mergarsal/FastCertRelPose/main/LICENSE.txt)


If you use this code for your research, please cite:
```
ADD OURS
```

Some parts of this repository are based on previous works: 
1. Matlab & python bindings from [TEASER++](https://arxiv.org/pdf/2001.07715.pdf)
        ```
                @article{yang2020teaser,
                  title={Teaser: Fast and certifiable point cloud registration},
                  author={Yang, Heng and Shi, Jingnan and Carlone, Luca},
                  journal={arXiv preprint arXiv:2001.07715},
                  year={2020}
                }
        ```
2. Scene generation from [opengv](https://github.com/laurentkneip/opengv.git) 



## Dependences 
* Eigen 
 ```
        sudo apt install libeigen3-dev
 ```

* Optimization (own fork)
 ```
        git clone https://github.com/mergarsal/Optimization.git
 ```
* GNCSO 
```
        git clone --recursive https://github.com/mergarsal/GNCSO.git 
```

* OpenGV
```
        git clone https://github.com/laurentkneip/opengv.git
```


## Build
```
git clone --recursive https://github.com/mergarsal/FastCertRelPose.git
cd GNCSO

mkdir build & cd build 

cmake .. 

make -jX

```

The compiled examples should be inside the `bin` directory. Run: 
```
        ./bin/example_essential_matrix
```
 


## Install 
In `build` folder: 
```
        sudo make install
```

We also provide the uninstall script: 
```
        sudo make uninstall
```



## How to use the library in your project

See the example in the folder `example_install` 
for the basic elements. 
       
1. In your CMakeLists.txt, add the dependences:
```
        find_package(gncso REQUIRED)
        find_package(opengv REQUIRED)        
        find_package(Essential REQUIRED)
```

2. For your executable, add the library in 
```
        target_link_libraries(XXX Essential opengv gncso)
```







