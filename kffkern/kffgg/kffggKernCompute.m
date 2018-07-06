function [K, Kbase, phid_X_hat] = kffggKernCompute(kern, X, X2)

% KFFGGKERNCOMPUTE Compute the KFFGG kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for
%	the KFF GG kernel given inputs associated with rows and
%	columns.
% RETURN K : the kernel matrix computed at the given points.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG X : the input matrix associated with the rows of the kernel.
% ARG X2 the input matrix associated with the columns of the kernel.
%
% FORMAT
% DESC computes the kernel matrix for the
%	KFF GG kernel given a design matrix of inputs.
%	 Returns:
% RETURN k : the kernel matrix computed at the given points.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG X : input data matrix in the form of a design matrix.
%
% SEEALSO : kffggKernParamInit, kernCompute, kernCreate, ggKernDiagCompute
%
% COPYRIGHT : Mauricio A. Alvarez, 2018

% KERN


if nargin < 3
    X2 = X;
end

S = kern.S;
Z = kern.Z;
Lambda = sqrt(kern.precisionU)*Z;
invPd = diag(1./kern.precisionG);
sqrtInvPd = sqrt(invPd);
LambdaSqrtInvPd    = Lambda*sqrtInvPd;
vLambdaSqrtInvPd   = sum(LambdaSqrtInvPd.^2,2);
mLambdaSqrtInvPdX  = repmat(vLambdaSqrtInvPd', size(X, 1), 1);
mLambdaSqrtInvPdX2 = repmat(vLambdaSqrtInvPd', size(X2, 1), 1);
phid_X  = exp(-0.5*mLambdaSqrtInvPdX  + 1j*X*Lambda');
phid_X2 = exp(-0.5*mLambdaSqrtInvPdX2 - 1j*X2*Lambda');
Kbase = phid_X*(phid_X2.');
K = (kern.sigma2Latent/S)*(kern.sensitivity^2)*real(Kbase);
if nargout > 2,
    phid_X_hat = sqrt(kern.sigma2Latent/S)*kern.sensitivity*phid_X;
end

