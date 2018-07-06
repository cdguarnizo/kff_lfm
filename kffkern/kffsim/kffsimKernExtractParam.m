function [params, names] = kffsimKernExtractParam(kern)

% KFFSIMKERNEXTRACTPARAM Extract params from the KFFSIM kernel structure.
% FORMAT
% DESC Extract parameters from the kernel Fourier features for sim into a 
% vector of parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. 
%
% FORMAT
% DESC Extract parameters from the kernel Fourier features for sim into a 
% vector of parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. 
% RETURN names : cell array of strings containing parameter names.
%
% SEEALSO kffsimKernParamInit, kffsimKernExpandParam
%
% COPYRIGHT : Mauricio A. Alvarez, 2017

% KERN

if nargout > 1
    [params, names] = simKernExtractParam(kern);
else
    params = simKernExtractParam(kern);
end