function [K, Kbase] = kffggXkffggKernCompute(ggKern1, ggKern2, X, X2)
% KFFGGXKFFGGKERNCOMPUTE Compute a cross kernel between two KFFGG kernels.
% FORMAT
% DESC computes cross kernel
%	terms between two KFF GG kernels for the multiple output kernel.
% RETURN K : block of values from kernel matrix.
% ARG ggkern1 : the kernel structure associated with the first KFF GG
% ARG ggkern2 : the kernel structure associated with the second KFF GG
%	   kernel.
% ARG x : inputs for which kernel is to be computed.
%
% DESC computes cross
%	kernel terms between two KFF GG kernels for the multiple output kernel.
% RETURN K :  block of values from kernel matrix.
% ARG ggkern1 : the kernel structure associated with the first KFF GG kernel.
% ARG ggkern2 : the kernel structure associated with the second KFF GG kernel.
% ARG x : row inputs for which kernel is to be computed.
% ARG x2 : column inputs for which kernel is to be computed.
%
% SEEALSO : multiKernParamInit, multiKernCompute, kffggKernParamInit
%
% COPYRIGHT : Mauricio A. Alvarez, 2018

% KERN

if nargin < 4
    X2 = X;
end

S = ggKern1.S;
Z = ggKern1.Z;
Lambda = sqrt(ggKern1.precisionU)*Z;

invPd = diag(1./ggKern1.precisionG);
sqrtInvPd = sqrt(invPd);
LambdaSqrtInvPd    = Lambda*sqrtInvPd;
vLambdaSqrtInvPd   = sum(LambdaSqrtInvPd.^2,2);
mLambdaSqrtInvPdX  = repmat(vLambdaSqrtInvPd', size(X, 1), 1);
phid_X  = (1/sqrt(S))*exp(-0.5*mLambdaSqrtInvPdX  + 1j*X*Lambda');

invPdp = diag(1./ggKern2.precisionG);
sqrtInvPdp = sqrt(invPdp);
LambdaSqrtInvPdp    = Lambda*sqrtInvPdp;
vLambdaSqrtInvPdp   = sum(LambdaSqrtInvPdp.^2,2);
mLambdaSqrtInvPdpX2  = repmat(vLambdaSqrtInvPdp', size(X2, 1), 1);
phidp_X2 = (1/sqrt(S))*exp(-0.5*mLambdaSqrtInvPdpX2 - 1j*X2*Lambda');

Kbase = phid_X*(phidp_X2.');
K = ggKern1.sigma2Latent*(ggKern1.sensitivity*ggKern2.sensitivity)*real(Kbase);


