function kern = kffgaussianKernExpandParam(kern, params)

% KFFGAUSSIANKERNEXPANDPARAM Create kernel structure from KFF Gaussian kernel's parameters.
% FORMAT
% DESC returns a KFF Gaussian kernel structure filled with the parameters in the given
%	vector. This is used as a helper function to enable parameters to be
%	optimised in, for example, the NETLAB optimisation functions.
% RETURN kern : kernel structure with the given parameters in the relevant
%	   locations.
% ARG kern : the kernel structure in which the parameters are to be
% ARG param : vector of parameters which are to be placed in the kernel
%	   structure.
%	
% SEEALSO : kffgaussianKernParamInit, kffgaussianKernExtractParam, kernExpandParam
%
% COPYRIGHT : Mauricio Alvarez, 2018

% KERN


kern.sigma2Latent = params(end);
kern.precisionU =  params(1:end-1)';