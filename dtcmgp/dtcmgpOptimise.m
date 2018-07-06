function [model, options, params] = dtcmgpOptimise(model, display, iters)

% MULTIGPOPTIMISE Optimise the inducing variable multigp based kernel.
% FORMAT
% DESC optimises the Gaussian
%	process  model for a given number of iterations.
% RETURN model : the optimised model.
% RETURN options : a vector containig the number of function and gradient
% RETURN params : the optimised parameter vector.
% evaluations. Also the number of iterations employed. 
% ARG model : the model to be optimised.
% ARG display : whether or not to display while optimisation proceeds,
%	   set to 2 for the most verbose and 0 for the least verbose.
% ARG iters : number of iterations for the optimisation.
%	
% SEEALSO : scg, conjgrad, multigpOptimiseCreate,
% multigpOptimiseGradient, multigpOptimiseObjective
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006
%
% MODIFIED : Mauricio A. Alvarez, 2008

% MULTIGP

if nargin < 3
  iters = 1000;
  if nargin < 2
    display = 1;
  end
end

params = dtcmgpExtractParam(model); 
options = optOptions;
if display,
    options(1) = 1;
     if length(params) <= 100 && display > 1
        options(9) = 1;
     end
end
options(14) = iters;

if isfield(model, 'optimiser'),
  optim = str2func(model.optimiser);
else
  optim = str2func('conjgrad');
end

if strcmp(func2str(optim), 'optimiMinimize')
  % Carl Rasmussen's minimize function 
  params = optim('dtcmgpObjectiveGradient', params, options, model);
elseif strcmp(func2str(optim), 'rprop')
    params = optim(params, 'dtcmgpObjectiveGradient', iters, model);
elseif strcmp(func2str(optim), 'scg2')
    [params, options] = optim('dtcmgpObjectiveGradient', params,  options, ...
                 'dtcmgpGradient', model);
else
  % NETLAB style optimization.
  [params, options] = optim('dtcmgpObjective', params,  options, ...
                 'dtcmgpGradient', model);
end
model = dtcmgpExpandParam(model, params);
model = dtcmgpComputeKernels(model);