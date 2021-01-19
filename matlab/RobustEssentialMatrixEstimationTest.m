close all; clear all; clc;


rng('shuffle');
% Test parameters:
% Default params.
noise=0.1;
FoV=150;
N = 200;  
% for t
min_parallax=0.5;
max_parallax=2.0;
% for the points
min_depth=1.0;
max_depth=8.0;

f = 800;  % focal length
% visualization


% weights vector
W = ones(N, 1);

% Params 
%  obsi = 0,
%  obsip,
%  weights,
%  max_res_tol_sq,
%  gnc_factor,
%  gnc_max_inner_iterations,
%  gnc_max_outer_iterations,
%  gnc_cost_threshold,
%  gnc_inliers_threshold,
%  gnc_min_nr_points,
%  init_method,
%  gnc_robust,
%  precon,
%  use_mult_certifiers,
%  verbose

max_res_tol_sq = 0.00001; 
gnc_factor = 1.10;
gnc_max_inner_iterations = 2; 
gnc_max_outer_iterations = 500;
gnc_cost_threshold = 0.000001;
gnc_inliers_threshold = 0.9; 
gnc_min_nr_points = 12;
gnc_robust = 5;

init_method = 0; % 8-pt
precon=2;        % dominant eigenvalues
use_mult_certifiers = 1; % check all the relaxations (\neq 6)
verbose=true; 

% 1. Generate data
[P1, P2, tgt, Rgt, indices_outliers] = create2D2DCorrespondencesNOutliers(N, noise, 80, FoV, min_parallax, max_parallax, ...
    min_depth, max_depth, false);

Egt = [[0 -tgt(3) tgt(2)]; [tgt(3) 0 -tgt(1)]; [-tgt(2) tgt(1) 0]]*Rgt;
Egt = Egt/Egt(3,3);




% Test the MEX function
[R, t, time_iterative, time_total, set_inliers, is_opt] = GNCEssentialMatrixMex(P1, P2, W', ...
                                                  max_res_tol_sq, gnc_factor, ...
                                                  gnc_max_inner_iterations, gnc_max_outer_iterations, ...
                                                  gnc_cost_threshold, gnc_inliers_threshold, ...
                                                  gnc_min_nr_points, init_method, gnc_robust, ...
                                                  precon, use_mult_certifiers, verbose);
% check error
dist_R(Rgt, R)



% Test the MEX wrapper function
[struct_output_wrap] = GNCEssentialMatrixEstimate(P1, P2, W', ...
                                                  'max_res_tol_sq', max_res_tol_sq, 'gnc_factor', gnc_factor, ...
                                                  'gnc_max_inner_iterations', gnc_max_inner_iterations, ...
                                                  'gnc_max_outer_iterations', gnc_max_outer_iterations, ...
                                                  'gnc_cost_threshold', gnc_cost_threshold, ...
                                                  'gnc_inliers_threshold', gnc_inliers_threshold, ...
                                                  'gnc_min_nr_points', gnc_min_nr_points, ...
                                                  'init_method', init_method, 'gnc_robust', gnc_robust, ...
                                                  'precon', precon, 'use_mult_certifiers', use_mult_certifiers, ...
                                                  'verbose', verbose);
dist_R(Rgt, struct_output_wrap.R)
