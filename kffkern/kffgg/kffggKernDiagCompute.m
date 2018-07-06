function K = kffggKernDiagCompute(kern, X)

% KFFGGKERNDIAGCOMPUTE Compute diagonal of KFF GG kernel.
% FORMAT
% DESC computes the diagonal of the kernel
%	matrix for the KFF Gaussian Gaussian kernel given a design matrix of
%	inputs.
% RETURN  K : a vector containing the diagonal of the kernel matrix computed
%	   at the given points.
% ARG kern : the kernel structure for which the matrix is computed.
% X : input data matrix in the form of a design matrix.
%	
% SEEALSO : kffggKernParamInit, kernDiagCompute, kernCreate, kffggKernCompute
%
% COPYRIGHT : Mauricio A. Alvarez, 2018

% KERN

S = kern.S;
Z = kern.Z;
Lambda = sqrt(kern.precisionU)*Z;
invPd = diag(1./kern.precisionG);
sqrtInvPd = sqrt(invPd);
LambdaSqrtInvPd    = Lambda*sqrtInvPd;
vLambdaSqrtInvPd   = sum(LambdaSqrtInvPd.^2,2);
phid_X  = exp(-0.5*vLambdaSqrtInvPd);
Kbase = phid_X.'*(phid_X);
K = (kern.sigma2Latent/S)*(kern.sensitivity^2)*real(Kbase)*ones(size(X,1),1);
