function k = kffgaussianKernDiagCompute(kern, x)

% KFFGAUSSIANKERNDIAGCOMPUTE Compute diagonal of KFF Gaussian kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the KFF Gaussian kernel
% given a design matrix of inputs.
% RETURN K : a vector containing the diagonal of the kernel matrix computed
% at the given points.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG X : input data matrix in the form of a design matrix.
%
% SEEALSO : kffgaussianKernParamInit, kernDiagCompute, kernCreate,
% gaussianKernCompute
%
% COPYRIGHT : Mauricio A. Alvarez, 2018

% KERN

k = kern.sigma2Latent*ones(size(x,1),1);
