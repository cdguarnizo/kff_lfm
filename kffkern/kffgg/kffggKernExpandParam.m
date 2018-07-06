function kern = kffggKernExpandParam(kern, params)

% KFFGGKERNEXPANDPARAM Create kernel structure from KFF GG kernel's parameters.
% FORMAT
% DESC returns a KFF Gaussian Gaussian
%	kernel structure filled with the parameters in the given vector.
%	This is used as a helper function to enable parameters to be
%	optimised in, for example, the NETLAB optimisation functions.
% RETURN kern : kernel structure with the given parameters in the relevant
%	   locations.
% ARG kern : the kernel structure in which the parameters are to be
%	   placed.
% ARG param : vector of parameters which are to be placed in the kernel
%	   structure.
%
% SEEALSO : kffggKernParamINit, kffggKernExtractParam, kernExpandParam
%
% COPYRIGHT : Mauricio A. Alvarez, 2018

% KERN

sizePU = size(kern.precisionU,1);
kern.precisionU = params(1:sizePU)';
sizePG = size(kern.precisionG,1);
kern.precisionG = params(sizePU+1:sizePU + sizePG)';
if ~(isfield(kern, 'isVarS') && kern.isVarS)    
    kern.sigma2Latent = params(end-1);
    kern.sensitivity = params(end);    
end


