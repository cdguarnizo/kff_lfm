function kern = kffggKernParamInit(kern, isArd)

% KFF GGKERNPARAMINIT KFF GG kernel parameter initialisation.
% FORMAT
% The KFF Gaussian Gaussian (gg) kernel corresponds to the covariance matrix of
% an output process y(x), which is generated through a convolution between a
% latent process with Gaussian covariance matrix u(s) and a Gaussian kernel
% K(x-s):
%	
%	y_n(x) =  sum_k int K_{nk}(x-s)u_k(s)ds
%	
% where K_{nk}(x-s) is a Gaussian kernel with precision matrix precision_y,
% and u_k(s) is an inducing function represented as a Gaussian process with
% inverse covariance precision_u.  With this assumptions, y_n(x) is also a
% Gaussian process with covariance provided by the Gaussian Gaussian kernel.
% Additionally, we use kernel Fourier features to approximate the kernel
% for u_k(s), that why we refer to this kernel as kff Gaussian Gaussian
%	
% The kernel is designed to interoperate with the multiple output block
% kernel so that u_k(s) can be inferred given several different
% instantiations of y_n(x).
%	
% Both precision_y and precision_u are considered as diagonal.
%
% DESC initialises the Gaussian Gaussian kernel structure with some default parameters.
% RETURN kern : the kernel structure with the default parameters placed in.
% ARG kern : the kernel structure which requires initialisation.
%	
% SEEALSO : kernCreate, kernParamInit, kffggKernCompute
%
% COPYRIGHT : Mauricio A. Alvarez, 2018

% KERN

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

kern.precisionU = 1;
if kern.isArd   
    kern.precisionG = ones(kern.inputDimension,1);
    kern.nParams = kern.inputDimension + 3;
else    
    kern.precisionG = 1;
    kern.nParams = 4;
end

kern.sigma2Latent = 1;
kern.sensitivity  = 1;
% The variances must be positive. As well as the sensitivity of the latent
% function.
kern.transforms.index = 1:(kern.nParams-1);
kern.transforms.type = optimiDefaultConstraint('positive');
kern.isVarS = false;
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
