function [struct_output] = EssentialMatrixEstimate(obsi, obsip, weights, varargin)



assert(size(obsi, 1) == 3, 'obsi must be a 3-by-N matrix.')
assert(size(obsip, 1) == 3, 'obsip must be a 3-by-N matrix.')
assert(size(weights, 1) == 1, 'weights must be a 1-by-N matrix.')


params = inputParser;
params.CaseSensitive = false;

addParameter(params, 'init_method', 0, ...
    @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=6);
    
addParameter(params, 'precon', 1, ...
    @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=4);
    
addParameter(params, 'use_mult_certifiers', 1, ...
    @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=6);
    
addParameter(params,'verbose', false, ...
    @(x) islogical(x));
    
parse(params, varargin{:});


[R, t, time_iterative, time_total, is_opt] = EssentialMatrixMex(obsi, obsip, weights, ...
                                     params.Results.init_method, ...
                                     params.Results.precon, ...
                                     params.Results.use_mult_certifiers, ...
                                     params.Results.verbose);

% struct_output = structEssentialEstimation;
struct_output = struct();
struct_output.R = R;
struct_output.t = t;
struct_output.time_iterative = time_iterative;
struct_output.time_total = time_total;
struct_output.is_opt = is_opt;

% class structEssentialEstimation
%   properties
%     R
%     t
%     time_iterative
%     time_total
%   end

%   end

  
end



