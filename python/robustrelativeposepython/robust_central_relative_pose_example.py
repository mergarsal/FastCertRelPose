#!/usr/bin/env python3.7

import numpy as np
import robustrelativeposepython
from numpy import linalg as LA



if __name__ == "__main__":
    

    # Generate random data points
    N = 20
    points3D = np.random.rand(3, N)

    # Apply arbitrary scale, translation and rotation
    translation = np.array([[1], [0.8], [-1]])
    rotation = np.array([[0.98370992, 0.17903344, -0.01618098],
                         [-0.04165862, 0.13947877, -0.98934839],
                         [-0.17486954, 0.9739059, 0.14466493]])
                         
    dst = np.matmul(rotation.transpose(), points3D - translation)
    
    # observations
    obs1 = points3D / LA.norm(points3D, axis=0)
    obs2 = dst / LA.norm(dst, axis=0)
    
    weights = np.ones(N)
    

    # set options 
    estimation_params = robustrelativeposepython.GNCEssentialClass.Options()
    estimation_params.estimation_verbose=0;
    estimation_params.gnc_robust = robustrelativeposepython.GNCEssentialClass.GNCRobustFunction.TUKEY;
    estimation_params.use_idx_relaxation=0;
    # create instance
    essential = robustrelativeposepython.GNCEssentialClass(obs1, obs2, weights,  
                estimation_params, np.ones([3,3])                )
    # solve
    results = essential.getResultGNC()
    # retrieve data
    
    # essential matrix
    print ("Essential matrix: ",     results.E_opt)
    # rotation matrix
    print ("Rotation matrix: ",     results.R_opt)
    print ("Gt rotation    : ",     rotation)
    # translation vector
    print ("Translation vector: ",     results.t_opt)
    print ("Gt Translation    : ",     translation / LA.norm(translation))
    # certifier status
    status = results.getStatus()
    print ("Status Estimation (0: optimal): ", status)
    # CPU time (total)
    print ("Elapsed time: ", results.elapsed_estimation_time)

