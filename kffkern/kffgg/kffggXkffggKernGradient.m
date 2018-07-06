function [g1, g2] = kffggXkffggKernGradient(ggKern1, ggKern2, X, X2, covGrad)

% KFFGGXKFFGGKERNGRADIENT Compute a cross gradient between two KFFGG kernels.
% FORMAT
% DESC computes cross gradient of parameters of a cross kernel between two
%	kff gg kernels for the multiple output kernel.
% RETURN g1 : gradient of the parameters of the first kernel, for ordering
%	   see kffggKernExtractParam.
% RETURN g2 : gradient of the parameters of the second kernel, for ordering
%	   see kffggKernExtractParam.
% ARG ggKern1 : the kernel structure associated with the first GG
%	   kernel.
% ARG ggKern2 : the kernel structure associated with the second GG
%	   kernel.
% ARG X : inputs for which kernel is to be computed.
% ARG covgrad : gradient of the objective function with respect to the
%	   elements of the cross kernel matrix.
%
% FORMAT
% DESC computes cross kernel terms between two KFF GG kernels for the multiple
%	output kernel.
% RETURN g1 : gradient of the parameters of the first kernel, for ordering
%	   see kffggKernExtractParam.
% RETURN g2 : gradient of the parameters of the second kernel, for ordering
%	   see kffggKernExtractParam.
% ARG ggKern1 : the kernel structure associated with the first GG
%	   kernel.
% ARG ggKern2 : the kernel structure associated with the second GG
%	   kernel.
% ARG X : row inputs for which kernel is to be computed.
% ARG X2 : column inputs for which kernel is to be computed.
% ARG covgrad : gradient of the objective function with respect to the
%	   elements of the cross kernel matrix.
%
% SEEALSO : multiKernParamInit, multiKernCompute, kffggKernParamInit,
% kffggKernExtractParam
%
% COPYRIGHT : Mauricio A. Alvarez, 2008

% KERN


if nargin < 5
    covGrad = X2;
    X2 = X;
end

S = ggKern1.S;
Z = ggKern1.Z;
Lambda = sqrt(ggKern1.precisionU)*Z;


Pd = diag(ggKern1.precisionG);
invPd = diag(1./ggKern1.precisionG);
sqrtInvPd = sqrt(invPd);
LambdaSqrtInvPd    = Lambda*sqrtInvPd;
vLambdaSqrtInvPd   = sum(LambdaSqrtInvPd.^2,2);
mLambdaSqrtInvPdX  = repmat(vLambdaSqrtInvPd', size(X, 1), 1);
phid_X  = exp(-0.5*mLambdaSqrtInvPdX  + 1j*X*Lambda');

Pdp = diag(ggKern2.precisionG);
invPdp = diag(1./ggKern2.precisionG);
sqrtInvPdp = sqrt(invPdp);
LambdaSqrtInvPdp    = Lambda*sqrtInvPdp;
vLambdaSqrtInvPdp   = sum(LambdaSqrtInvPdp.^2,2);
mLambdaSqrtInvPdpX2  = repmat(vLambdaSqrtInvPdp', size(X2, 1), 1);
phidp_X2 = exp(-0.5*mLambdaSqrtInvPdpX2 - 1j*X2*Lambda');
Kbase = phid_X*(phidp_X2.');
Lambda2 = Lambda.^2;

% Gradient wrt Pd and Pdp
matGradPd = zeros(ggKern1.inputDimension, 1);
matGradPdp = zeros(ggKern2.inputDimension, 1);
for i=1:ggKern1.inputDimension
    phid_X_dot_Lambda = phid_X.*repmat((Lambda2(:, i))', size(X, 1), 1);
    gradPHIXPHIX2Pd = (0.5/(Pd(i,i)^2))*(phid_X_dot_Lambda*(phidp_X2.'));
    matGradPd(i) = sum(sum(covGrad.*real(gradPHIXPHIX2Pd)));   
    phid_X2_dot_Lambda = phidp_X2.*repmat((Lambda2(:, i))', size(X2, 1), 1);
    gradPHIXPHIX2Pdp = (0.5/(Pdp(i,i)^2))*(phid_X*(phid_X2_dot_Lambda.'));
    matGradPdp(i) = sum(sum(covGrad.*real(gradPHIXPHIX2Pdp)));   
end
matGradPd  = (ggKern1.sigma2Latent/S)*(ggKern1.sensitivity*ggKern2.sensitivity)*matGradPd;
matGradPdp = (ggKern1.sigma2Latent/S)*(ggKern1.sensitivity*ggKern2.sensitivity)*matGradPdp;

% Gradient wrt precisionU

ZSqrtInvPd    = Z*sqrtInvPd;
vZSqrtInvPd   = sum(ZSqrtInvPd.^2,2);
mZSqrtInvPdX  = repmat(vZSqrtInvPd', size(X, 1), 1);
ZSqrtInvPdp    = Z*sqrtInvPdp;
vZSqrtInvPdp   = sum(ZSqrtInvPdp.^2,2);
mZSqrtInvPdpX2 = repmat(vZSqrtInvPdp', size(X2, 1), 1);
CX  = X*Z';
CX2 = X2*Z';

gradPHIXPHIX2 = -0.5*((phid_X.*mZSqrtInvPdX)*phidp_X2.') ...
                +(0.5*1i)*(ggKern1.precisionU^(-1/2))*((phid_X.*CX)*phidp_X2.') ...  
                -0.5*(phid_X*(phidp_X2.*mZSqrtInvPdpX2).') ...
                - (0.5*1i)*(ggKern1.precisionU^(-1/2))*(phid_X*(phidp_X2.*CX2).');

matGradPrecU = real((ggKern1.sigma2Latent/S)*(ggKern1.sensitivity*ggKern2.sensitivity)*...
    sum(sum(covGrad.*gradPHIXPHIX2)));

% Gradients sigma2 and sensitivity

gradconst = real(sum(sum(covGrad.*Kbase)));
gradSigma2Latent = (ggKern1.sensitivity*ggKern2.sensitivity/S)*gradconst;
gradSens1 = (ggKern1.sigma2Latent*ggKern2.sensitivity/S)*gradconst;
gradSens2 = (ggKern1.sigma2Latent*ggKern1.sensitivity/S)*gradconst;

%%%%%

g1 = [matGradPrecU  matGradPd(:)'   gradSigma2Latent gradSens1];
g2 = [0             matGradPdp(:)'  0                gradSens2];
 
