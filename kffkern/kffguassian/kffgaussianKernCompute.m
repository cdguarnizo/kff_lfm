function [K, Kbase, phi_X, phi_X2]  = kffgaussianKernCompute(kern, X, X2)

% KFFGAUSSIANKERNCOMPUTE Compute the Gaussian kernel with Fourier Features.
% FORMAT
% DESC computes the kernel parameters for the KFF Gaussian kernel given
% inputs associated with rows and columns.
% RETURN K : the kernel matrix computed at the given points.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG X : the input matrix associated with the rows of the kernel.
% ARG X2 : the input matrix associated with the columns of the kernel.
%
% FORMAT
% DESC computes the kernel matrix for the KFF Gaussian kernel given a design matrix of inputs.
% RETURN K : the kernel matrix computed at the given points.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG X : input data matrix in the form of a design matrix.
%
% SEEALSO : kffgaussianKernParamInit, kernCompute, kernCreate, kffgaussianKernDiagCompute
%
% COPYRIGHT : Mauricio Alvarez, 2018
%

% KERN

if nargin < 3
    X2 = X;
end
S = kern.S;
Z = kern.Z;
Lambda = sqrt(kern.precisionU)*Z;

phi_X  = exp(1i*X*Lambda');
phi_X2 = exp(-1i*X2*Lambda');

Kbase = phi_X*(phi_X2.');
K = (kern.sigma2Latent/S)*real(Kbase);







