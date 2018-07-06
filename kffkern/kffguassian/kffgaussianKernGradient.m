function g = kffgaussianKernGradient(kern, X, varargin)

% KFFGAUSSIANKERNGRADIENT Gradient of KFF Gaussian kernel's parameters.
% FORMAT
% DESC computes the gradient of
%	functions with respect to the KFF Gaussian kernel's
%	parameters. As well as the kernel structure and the input positions,
%	the user provides a matrix PARTIAL which gives the partial
%	derivatives of the function with respect to the relevant elements of
%	the kernel matrix.
% RETURN g:  gradients of the function of interest with respect to the
%	   kernel parameters. The ordering of the vector should match that
%	   provided by the function kernExtractParam.
% ARG kern : the kernel structure for which the gradients are being
%	   computed.
% ARG X : the input locations for which the gradients are being
%	   computed.
% ARG covGrad : matrix of partial derivatives of the function of
%	   interest with respect to the kernel matrix. The argument takes the
%	   form of a square matrix of dimension  numData, where numData is
%	   the number of rows in X.
%
% FORMAT
% DESC  computes the derivatives
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
% ARG covGrad : matrix of partial derivatives of the function of
%	   interest with respect to the kernel matrix. The matrix should have
%	   the same number of rows as X1 and the same number of columns as X2
%	   has rows.
%
% SEEALSO : kffgaussianKernParamInit, kernGradient,
%
% COPYRIGHT : Mauricio A. Alvarez, 2018

% KERN

if nargin < 4
    X2 = X;
    covGrad = varargin{1};
else
    X2 = varargin{1};
    covGrad = varargin{2};
end

S = kern.S;
Z = kern.Z;
Lambda = sqrt(kern.precisionU)*Z;

phi_X  = exp(1i*X*Lambda');
phi_X2 = exp(-1i*X2*Lambda');

Kbase = real(phi_X*(phi_X2.'));

CX  = X*Z';
CX2 = X2*Z';

gradPHIXPHIX2 = (0.5*1i)*(kern.precisionU^(-1/2))*((phi_X.*CX)*phi_X2.') ...    
    - (0.5*1i)*(kern.precisionU^(-1/2))*(phi_X*(phi_X2.*CX2).');

matGrad = real((kern.sigma2Latent/S)*sum(sum(covGrad.*gradPHIXPHIX2)));

g = [matGrad(:)' sum(sum(covGrad.*Kbase))/S];



