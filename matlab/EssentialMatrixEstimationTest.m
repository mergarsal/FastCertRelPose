close all; clear all; clc;


rng('shuffle');
% Test parameters:
% Default params.
noise=1.5;
FoV=100;
N = 12;
N = N;
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

% 1. Generate data
[P1, P2, tgt, Rgt, indices_outliers] = create2D2DCorrespondencesNOutliers(N, noise, 0, FoV, min_parallax, max_parallax, ...
    min_depth, max_depth, false);

Egt = [[0 -tgt(3) tgt(2)]; [tgt(3) 0 -tgt(1)]; [-tgt(2) tgt(1) 0]]*Rgt;
Egt = Egt/Egt(3,3);


%% Check original fcns:

% 
tic
[E_mani, is_opt_mani, optval_mani, dualval, dual_gap] = estimateCertifiedEssentialMatrix(P1, P2, W);
toc


tic
[E_mani, is_opt_mani, optval_mani, dualval, dual_gap] = estimateCertNaiveManifold(P1, P2, W);
toc;



% Test the MEX function
% enum class INPUT_PARAMS : int
%{
%  obsi = 0,
%  obsip,
%  weights,
%  init_method,
%  precon,
%  use_mult_certifiers,
%  verbose
%};
%
%enum class OUTPUT_PARAMS : int
%{
%  R_est = 0,
%  t_est,
%  time_total,
%  time_iterative,
%  is_opt
%};

init_method = 0; % 8-pt
precon=2;        % dominant eigenvalues
use_mult_certifiers = 1; % check all the relaxations (\neq 6)
verbose=true; 

[R, t, time_iterative, time_total, is_opt] = EssentialMatrixMex(P1, P2, W', ...
        init_method, precon, use_mult_certifiers, verbose);
% check error
dist_R(Rgt, R)



% Test the MEX wrapper function
[struct_output_wrap] = EssentialMatrixEstimate(P1, P2, W', ...
        'init_method', init_method, 'precon', precon, 'use_mult_certifiers', use_mult_certifiers, 'verbose', verbose);
dist_R(Rgt, struct_output_wrap.R)
