function [g1, g2] = kffggXkffgaussianKernGradient(ggKern, gaussianKern, X, X2, covGrad)
% KFFGGXKFFGAUSSIANKERNGRADIENT Compute gradient between the KFF GG and KFF GAUSSIAN kernels.
% FORMAT
% DESC computes the
%	gradient of an objective function with respect to cross kernel terms
%	between KFF GG and KFF GAUSSIAN kernels for the multiple output kernel.
% RETURN g1 : gradient of objective function with respect to kernel
%	   parameters of KFF GG kernel.
% RETURN g2 : gradient of objective function with respect to kernel
%	   parameters of KFF GAUSSIAN kernel.
% ARG ggkern : the kernel structure associated with the GG kernel.
% ARG gaussianKern :  the kernel structure associated with the GAUSSIAN kernel.
% ARG x : inputs for which kernel is to be computed.
%
% FORMAT
% DESC  computes
%	the gradient of an objective function with respect to cross kernel
%	terms between KFF GG and KFF GAUSSIAN kernels for the multiple output kernel.
% RETURN g1 : gradient of objective function with respect to kernel
%	   parameters of KFF GG kernel.
% RETURN g2 : gradient of objective function with respect to kernel
%	   parameters of KFF GAUSSIAN kernel.
% ARG ggKern : the kernel structure associated with the KFF GG kernel.
% ARG gaussianKern : the kernel structure associated with the GAUSSIAN kernel.
% ARG x1 : row inputs for which kernel is to be computed.
% ARG x2 : column inputs for which kernel is to be computed.
%
% SEEALSO : multiKernParamInit, multiKernCompute, kffggKernParamInit,
% gaussianKernParamInit
%
% COPYRIGHT : Mauricio A. Alvarez and Neil D. Lawrence, 2008
%
% MODIFICATIONS : Mauricio A. Alvarez, 2009, 2013
 
% KERN
 
if nargin < 5
    covGrad = X2;
    X2 = X;
end

S = ggKern.S;
Z = ggKern.Z;
Lambda = sqrt(ggKern.precisionU)*Z;

Pd = diag(ggKern.precisionG);
invPd = diag(1./ggKern.precisionG);
sqrtInvPd = sqrt(invPd);
LambdaSqrtInvPd    = Lambda*sqrtInvPd;
vLambdaSqrtInvPd   = sum(LambdaSqrtInvPd.^2,2);
mLambdaSqrtInvPdX  = repmat(vLambdaSqrtInvPd', size(X, 1), 1);
phid_X  = exp(-0.5*mLambdaSqrtInvPdX  + 1j*X*Lambda');
xi_X2   = exp(-1j*X2*Lambda');
Kbase = phid_X*(xi_X2.');

Lambda2 = Lambda.^2;

% Gradient wrt Pd 
matGradPd = zeros(ggKern.inputDimension, 1);
for i=1:ggKern.inputDimension
    phid_X_dot_Lambda = phid_X.*repmat((Lambda2(:, i))', size(X, 1), 1);
    gradPHIXPHIX2Pd = (0.5/(Pd(i,i)^2))*(phid_X_dot_Lambda*(xi_X2.'));
    matGradPd(i) = sum(sum(covGrad.*real(gradPHIXPHIX2Pd)));   
end
matGradPd  = (ggKern.sigma2Latent/S)*(ggKern.sensitivity)*matGradPd;

% Gradient wrt precisionU 

ZSqrtInvPd    = Z*sqrtInvPd;
vZSqrtInvPd   = sum(ZSqrtInvPd.^2,2);
mZSqrtInvPdX  = repmat(vZSqrtInvPd', size(X, 1), 1);
CX  = X*Z';
CX2 = X2*Z';

gradPHIXPHIX2 = -0.5*((phid_X.*mZSqrtInvPdX)*xi_X2.') ...
                +(0.5*1i)*(ggKern.precisionU^(-1/2))*((phid_X.*CX)*xi_X2.') ...                  
                - (0.5*1i)*(ggKern.precisionU^(-1/2))*(phid_X*(xi_X2.*CX2).');

matGradPrecU = real((ggKern.sigma2Latent/S)*(ggKern.sensitivity)*...
    sum(sum(covGrad.*gradPHIXPHIX2)));

% Gradients sigma2 and sensitivity

gradconst = real(sum(sum(covGrad.*Kbase)));
gradSigma2Latent = (ggKern.sensitivity/S)*gradconst;
gradSens1 = (gaussianKern.sigma2Latent/S)*gradconst;

g1 = [matGradPrecU matGradPd(:)' gradSigma2Latent gradSens1];
g2 = [0 0];
