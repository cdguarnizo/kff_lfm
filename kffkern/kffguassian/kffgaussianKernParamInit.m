function kern = kffgaussianKernParamInit(kern)

% KFFGAUSSIANKERNPARAMINIT Gaussian kernel parameter initialisation with
% Fourier features
% The Gaussian kernel used here follows the shape of a gaussian
%	distribution computed using Fourier Features
%
%	k(x_i, x_j) =  sigma2*(1/S)*sum_{\forall s} exp(- 0.5*lambda_s(x_i - x_j))
%
%	In the above equation, lambda_s has been sampled from a Gaussian with precisionU 
%   and sigma2 is a variance factor. 
%
%  FORMAT
% DESC  initialises the KFF Gaussian
%	kernel structure with some default parameters.
% RETURN kern : the kernel structure with the default parameters placed in.
% ARG kern : the kernel structure which requires initialisation.
%
%
% SEEALSO : kernCreate, kernParamInit
%
% COPYRIGHT : Mauricio A. Alvarez, 2018
%

% KERN

if kern.inputDimension == 0
    kern.inputDimension = 1;
end

if isfield(kern, 'options') && isfield(kern.options, 'isArd')
    kern.isArd = kern.options.isArd;
else
    switch nargin
        case 1
            kern.isArd = false;
        case 2
            kern.isArd = isArd;
    end
end

kern.sigma2Latent = 1;
if kern.isArd
    % To be developed
else
    kern.precisionU   = 1;
    kern.nParams = 2;
end
kern.isVarS = false;
% Constrains parameters positive for optimisation.
% The variances of P need to be positive and we constrain the sensitivity
% to be positive as well
kern.transforms.index =1:kern.nParams;
kern.transforms.type = optimiDefaultConstraint('positive');
kern.isStationary = true;

if isfield(kern, 'options') && isfield(kern.options, 'S') 
    kern.S = kern.options.S;    
else
    kern.S = 10;
end
if isfield(kern, 'options') && isfield(kern.options, 'Z') 
    kern.Z = kern.options.Z;    
else
    kern.Z = randn(kern.S, kern.inputDimension); % Create the basis functions in the kernel Fourier feature way
end
