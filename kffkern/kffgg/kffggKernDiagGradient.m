function g = kffggKernDiagGradient(kern, X, covGrad)

% KFFGGKERNDIAGGRADIENT Compute gradient for the diagonal of KFF GG kernel.
% FORMAT
% DESC computes the gradient for the diagonal of the kernel
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
Pd = diag(kern.precisionG);
invPd = diag(1./kern.precisionG);
sqrtInvPd = sqrt(invPd);
LambdaSqrtInvPd    = Lambda*sqrtInvPd;
vLambdaSqrtInvPd   = sum(LambdaSqrtInvPd.^2,2);
phid_X  = exp(-0.5*vLambdaSqrtInvPd);
Kbase = real(phid_X.'*(phid_X))*ones(size(X,1),1);

% Gradient wrt Pd
matGradPd = zeros(kern.inputDimension,1);
for i=1:kern.inputDimension
    phid_X_dot_Lambda = phid_X.*Lambda(:, i);    
    gradPHIXPHIX = (1/(Pd(i,i)^2))*(phid_X_dot_Lambda.'*phid_X_dot_Lambda);
    matGradPd(i) = sum(covGrad.*real(gradPHIXPHIX*ones(size(X,1),1)));   
end
matGradPd = (kern.sigma2Latent/S)*(kern.sensitivity^2)*matGradPd;

%Gradient wrt precisionU
ZSqrtInvPd    = Z*sqrtInvPd;
vZSqrtInvPdX   = sum(ZSqrtInvPd.^2,2);

gradPHIXPHIX = -0.5*((phid_X.*vZSqrtInvPdX).'*phid_X) ...                 
                -0.5*(phid_X.'*(phid_X.*vZSqrtInvPdX));
            
matGradPrecU = real((kern.sigma2Latent/S)*(kern.sensitivity^2)*...
    sum(covGrad.*(gradPHIXPHIX*ones(size(X,1),1))));

% Gradients sigma2 and sensitivity

gradconst = real((sum(covGrad.*Kbase)));
gradSigma2Latent = (kern.sensitivity^2/S)*gradconst;
gradSensitivity = (2*kern.sigma2Latent*kern.sensitivity/S)*gradconst;

g = [matGradPrecU(:)' matGradPd(:)' gradSigma2Latent gradSensitivity];

