function [K, Kbase] = kffggXkffgaussianKernCompute(ggKern, gaussianKern, X, X2)

% KFFGGXKFFGAUSSIANKERNCOMPUTE Compute a cross kernel between the KFF GG and KFF GAUSSIAN kernels.
% FORMAT
% DESC computes cross kernel
%	terms between KFF GG and KFF GAUSSIAN kernels for the multiple output kernel.
% RETURN K :  block of values from kernel matrix.
% ARG ggKern : the kernel structure associated with the KFF GG kernel.
% ARG gaussianKern : the kernel structure associated with the KFF GAUSSIAN kernel.
% ARG X :  inputs for which kernel is to be computed.
%
% FORMAT
% DESC computes cross
%	kernel terms between KFF GG and KFF GAUSSIAN kernels for the multiple output
%	kernel.
% RETURN K : block of values from kernel matrix.
% ARG ggKern :  the kernel structure associated with the KFF GG kernel.
% ARG gaussianKern the kernel structure associated with the KFF GAUSSIAN kernel.
% ARG X : row inputs for which kernel is to be computed.
% ARG X2 : column inputs for which kernel is to be computed.
%
% SEEALSO : multiKernParamInit, multiKernCompute, kffggKernParamInit, kffgaussianKernParamInit
%
%
% COPYRIGHT : Mauricio A. Alvarez, 2018

% KERN

if nargin < 4
    X2 = X;
end

S = ggKern.S;
Z = ggKern.Z;
Lambda = sqrt(ggKern.precisionU)*Z;

invPd = diag(1./ggKern.precisionG);
sqrtInvPd = sqrt(invPd);
LambdaSqrtInvPd    = Lambda*sqrtInvPd;
vLambdaSqrtInvPd   = sum(LambdaSqrtInvPd.^2,2);
mLambdaSqrtInvPdX  = repmat(vLambdaSqrtInvPd', size(X, 1), 1);
phid_X  = exp(-0.5*mLambdaSqrtInvPdX  + 1j*X*Lambda');
xi_X2   = exp(-1j*X2*Lambda');
Kbase = phid_X*(xi_X2.');
K = (gaussianKern.sigma2Latent/S)*ggKern.sensitivity*real(Kbase);




