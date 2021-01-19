function [struct_output] = GNCEssentialMatrixEstimate(obsi, obsip, weights, varargin)
%% Add info here


assert(size(obsi, 1) == 3, 'obsi must be a 3-by-N matrix.')
assert(size(obsip, 1) == 3, 'obsip must be a 3-by-N matrix.')
assert(size(weights, 1) == 1, 'weights must be a 1-by-N matrix.')

% enum class INPUT_PARAMS : int
% {
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
%};


params = inputParser;
params.CaseSensitive = false;

addParameter(params, 'max_res_tol_sq', 0.01^2, ...
    @(x) isnumeric(x) && isscalar(x) && x>=0);
    
addParameter(params, 'gnc_factor', 1.10, ...
    @(x) isnumeric(x) && isscalar(x) && x>1);
    
addParameter(params, 'gnc_max_inner_iterations', 2, ...
    @(x) isnumeric(x) && isscalar(x) && x>=1);
    
addParameter(params, 'gnc_max_outer_iterations', 500, ...
    @(x) isnumeric(x) && isscalar(x) && x>=1);
    
addParameter(params, 'gnc_cost_threshold', 0.0000001, ...
    @(x) isnumeric(x) && isscalar(x) && x>=0);
    
addParameter(params, 'gnc_inliers_threshold', 0.9, ...
    @(x) isnumeric(x) && isscalar(x) && x>0 && x<1);
    
addParameter(params, 'gnc_min_nr_points', 12, ...
    @(x) isnumeric(x) && isscalar(x) && x>0 );
    
addParameter(params, 'init_method', 0, ...
    @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=6);
    
addParameter(params, 'gnc_robust', 2, ...
    @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=6);
    
addParameter(params, 'precon', 1, ...
    @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=4);
    
addParameter(params, 'use_mult_certifiers', 1, ...
    @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=6);
    
addParameter(params,'verbose', false, ...
    @(x) islogical(x));
    
parse(params, varargin{:});


[R, t, time_iterative, time_total, set_inliers, is_opt] = GNCEssentialMatrixMex(obsi, obsip, weights, ...
                                                  params.Results.max_res_tol_sq, params.Results.gnc_factor, ...
                                                  params.Results.gnc_max_inner_iterations, ...
                                                  params.Results.gnc_max_outer_iterations, params.Results.gnc_cost_threshold, ...
                                                  params.Results.gnc_inliers_threshold, params.Results.gnc_min_nr_points, ...
                                                  params.Results.init_method, params.Results.gnc_robust, ...
                                                  params.Results.precon, ...
                                                  params.Results.use_mult_certifiers, params.Results.verbose);

struct_output = struct();
struct_output.R = R;
struct_output.t = t;
struct_output.time_iterative = time_iterative;
struct_output.time_total = time_total;
struct_output.set_inliers = set_inliers;
struct_output.is_opt = is_opt;

end
