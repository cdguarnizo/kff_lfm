function kern = kffsimKernExpandParam(kern, params)

% KFFSIMKERNEXPANDPARAM Kernel structure from KFFSIM kernel's parameters.
% FORMAT
% DESC returns a kernel Fourier feature SIM structure filled with the
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
% SEEALSO : simKernParamInit, simKernExtractParam, kernExpandParam
%
% COPYRIGHT : Mauricio A Alvarez, 2017

% KERN

kern = simKernExpandParam(kern, params);