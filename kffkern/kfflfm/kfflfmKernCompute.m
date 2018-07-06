function k = kfflfmKernCompute(kern, t, t2)

% KFFLFMKERNCOMPUTE Compute the KFFLFM kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for the KFF latent force model
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG t1 : the input matrix associated with the rows of the kernel.
% ARG t2 : the input matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel matrix for the KFF latent force model
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG t : input data matrix in the form of a design matrix.
% RETURN k : the kernel matrix computed at the given points.
%
% SEEALSO : lfmKernParamInit, kernCompute, kernCreate, lfmKernDiagCompute
%
% COPYRIGHT : Mauricio A. Alvarez, 2017

% KERN


if nargin < 3
  t2 = t;
end
if size(t, 2) > 1 || size(t2, 2) > 1
  error('Input can only have one column');
end

k = kfflfmXkfflfmKernCompute(kern, kern, t, t2);

if nargin < 3
  k = k + k';
  k = k*0.5;
end

k = real(k); % introduced MA 2008
