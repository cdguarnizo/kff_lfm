function kern = kfflfmKernExpandParam(kern, params)

% KFFLFMKERNEXPANDPARAM Kernel structure from KFFLFM kernel's parameters.
% FORMAT
% DESC returns a kernel Fourier feature LFM structure filled with the
% parameters in the given vector. This is used as a helper function to
% enable parameters to be optimised in, for example, the NETLAB
% optimisation functions.
% ARG kern : the kernel structure in which the parameters are to be
% placed.
% ARG param : vector of parameters which are to be placed in the
% kernel structure.
% RETURN kern : kernel structure with the given parameters in the
% relevant locations.
%
% SEEALSO : lfmKernParamInit, lfmKernExtractParam, kernExpandParam
%
% COPYRIGHT : Mauricio A Alvarez, 2017

% KERN

kern = lfmKernExpandParam(kern, params);