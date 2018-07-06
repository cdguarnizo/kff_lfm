function g = kffggKernGradient(kern, X, varargin)

% KFFGGKERNGRADIENT Gradient of KFF GG kernel's parameters.
% FORMAT
% DESC computes the gradient of
%	functions with respect to the KFF Gaussian Gaussian kernel's
%	parameters. As well as the kernel structure and the input positions,
%	the user provides a matrix PARTIAL which gives the partial
%	derivatives of the function with respect to the relevant elements of
%	the kernel matrix.
% RETURN g : gradients of the function of interest with respect to the
%	   kernel parameters. The ordering of the vector should match that
%	   provided by the function kernExtractParam.
% ARG kern : the kernel structure for which the gradients are being
%	   computed.
% ARG X : the input locations for which the gradients are being
%	   computed.
% ARG partial : matrix of partial derivatives of the function of
%	   interest with respect to the kernel matrix. The argument takes the
%	   form of a square matrix of dimension  numData, where numData is
%	   the number of rows in X.
%
% FORMAT
% DESC computes the derivatives
%	as above, but input locations are now provided in two matrices
%	associated with rows and columns of the kernel matrix.
% RETURN g : gradients of the function of interest with respect to the
%	   kernel parameters.
% ARG kern : the kernel structure for which the gradients are being
%	   computed.
% ARG X1 : the input locations associated with the rows of the kernel
%	   matrix.
% ARG X2 : the input locations associated with the columns of the kernel
%	   matrix.
% ARG partial : matrix of partial derivatives of the function of
%	   interest with respect to the kernel matrix. The matrix should have
%	   the same number of rows as X1 and the same number of columns as X2
%	   has rows.
%
% SEEALSO : kffggKernParamInit, kffkernGradient, ggKernDiagGradient, kernGradX
%
% COPYRIGHT : Mauricio A. Alvarez, 2018

% KERN

if length(varargin)<2
    X2 = X;
else
    X2 = varargin{1};
end

covGrad = varargin{end};

S = kern.S;
Z = kern.Z;
Lambda = sqrt(kern.precisionU)*Z;
Pd  = diag(kern.precisionG);
invPd = diag(1./kern.precisionG);
sqrtInvPd = sqrt(invPd);
LambdaSqrtInvPd    = Lambda*sqrtInvPd;
vLambdaSqrtInvPd   = sum(LambdaSqrtInvPd.^2,2);
mLambdaSqrtInvPdX  = repmat(vLambdaSqrtInvPd', size(X, 1), 1);
mLambdaSqrtInvPdX2 = repmat(vLambdaSqrtInvPd', size(X2, 1), 1);
phid_X  = exp(-0.5*mLambdaSqrtInvPdX  + 1j*X*Lambda');
phid_X2 = exp(-0.5*mLambdaSqrtInvPdX2 - 1j*X2*Lambda');
Kbase = phid_X*(phid_X2.');

% Gradient wrt Pd
matGradPd = zeros(kern.inputDimension,1);
for i=1:kern.inputDimension
    phid_X_dot_Lambda = phid_X.*repmat((Lambda(:, i))', size(X, 1), 1);
    phid_X2_dot_Lambda = phid_X2.*repmat((Lambda(:, i))', size(X2, 1), 1);
    gradPHIXPHIX2 = (1/(Pd(i,i)^2))*(phid_X_dot_Lambda*(phid_X2_dot_Lambda.'));
    matGradPd(i) = sum(sum(covGrad.*real(gradPHIXPHIX2)));   
end
matGradPd = (kern.sigma2Latent/S)*(kern.sensitivity^2)*matGradPd;

%Gradient wrt precisionU

ZSqrtInvPd    = Z*sqrtInvPd;
vZSqrtInvPd   = sum(ZSqrtInvPd.^2,2);
mZSqrtInvPdX  = repmat(vZSqrtInvPd', size(X, 1), 1);
mZSqrtInvPdX2 = repmat(vZSqrtInvPd', size(X2, 1), 1);
CX  = X*Z';
CX2 = X2*Z';

gradPHIXPHIX2 = -0.5*((phid_X.*mZSqrtInvPdX)*phid_X2.') ...
                +(0.5*1i)*(kern.precisionU^(-1/2))*((phid_X.*CX)*phid_X2.') ...  
                -0.5*(phid_X*(phid_X2.*mZSqrtInvPdX2).') ...
                - (0.5*1i)*(kern.precisionU^(-1/2))*(phid_X*(phid_X2.*CX2).');

matGradPrecU = real((kern.sigma2Latent/S)*(kern.sensitivity^2)*...
    sum(sum(covGrad.*gradPHIXPHIX2)));

% Gradients sigma2 and sensitivity

gradconst = real(sum(sum(covGrad.*Kbase)));
gradSigma2Latent = (kern.sensitivity^2/S)*gradconst;
gradSensitivity = (2*kern.sigma2Latent*kern.sensitivity/S)*gradconst;

g = [matGradPrecU(:)' matGradPd(:)' gradSigma2Latent gradSensitivity];
