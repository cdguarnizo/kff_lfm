function g = kffsimKernGradient(kern, t, varargin)

% KFFSIMKERNGRADIENT Gradient of KFFSIM kernel's parameters.
% FORMAT
% DESC computes the kernel parameters for the kff single input motif
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG t1 : the input matrix associated with the rows of the kernel.
% ARG t2 : the input matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
% RETURN sk : unscaled kernel matrix (i.e. only 0.5 times the sum of h's
% part).
%
% FORMAT
% DESC computes the kernel matrix for the kff single input motif
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG t : input data matrix in the form of a design matrix.
% RETURN k : the kernel matrix computed at the given points.
% RETURN sk : unscaled kernel matrix (i.e. only 0.5 times the sum of h's
% part).
%
% SEEALSO : kffsimKernParamInit, kernCompute, kernCreate, simKernDiagCompute
%
% COPYRIGHT :  Mauricio Alvarez, 2017


if length(varargin)<2
  t2 = t;
else
  t2 = varargin{1};
end

[g1, g2] = kffsimXkffsimKernGradient(kern, kern, t, t2, varargin{end});

g = g1 + g2;

if isfield(kern, 'gaussianInitial') && kern.gaussianInitial, %TODO: Fix this code for the conditional
  dim1 = size(t, 1);
  dim2 = size(t2, 1);
  t1Mat = t(:, ones(1, dim2));
  t2Mat = t2(:, ones(1, dim1))';

  % Gradient with respect to the decay
  g(1) = g(1) + sum(sum(-kern.initialVariance*(t1Mat + t2Mat).*exp(-kern.decay*(t1Mat + t2Mat)) .* varargin{end}));
  g(end+1) = sum(sum(exp(-kern.decay*(t1Mat + t2Mat)) .* varargin{end}));
end
