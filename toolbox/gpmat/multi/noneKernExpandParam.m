function kern = noneKernExpandParam(kern, params)

% NONEKERNEXPANDPARAM Create kernel structure from NONE kernel's parameters.
% FORMAT
% DESC returns a dummy kernel function kernel structure filled with the
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
% SEEALSO : noneKernParamInit, noneKernExtractParam, kernExpandParam
%
% COPYRIGHT : Neil D. Lawrence, 2008

% KERN

kern = kern;