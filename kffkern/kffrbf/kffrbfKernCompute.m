function K = kffrbfKernCompute(lfmKern1, t1, t2)

% LFMXLFMKERNCOMPUTE Compute a cross kernel between two KFF LFM kernels.
% FORMAT
% DESC computes cross kernel terms between two LFM kernels for
% the multiple output kernel.
% ARG lfmKern1 : the kernel structure associated with the first LFM
% kernel.
% ARG lfmKern2 : the kernel structure associated with the second LFM
% kernel.
% ARG t : inputs for which kernel is to be computed.
% RETURN K : block of values from kernel matrix.
%
% FORMAT
% DESC computes cross kernel terms between two KFF LFM kernels for
% the multiple output kernel.
% ARG lfmKern1 : the kernel structure associated with the first LFM
% kernel.
% ARG lfmKern2 : the kernel structure associated with the second LFM
% kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% RETURN K : block of values from kernel matrix.
%
% SEEALSO : lfmKernParamInit, lfmKernCompute, lfmKernParamInit
%
% COPYRIGHT : Mauricio Alvarez, 2017

% KERN

if nargin < 3,
    t2 = t1;
end
if size(t1, 2) > 1 || size(t2, 2) > 1,
    error('Input can only have one column');
end

% using complex numbers
S = lfmKern1.S; 
Z = lfmKern1.Z;

lambda_s = sqrt(lfmKern1.inverseWidth)*Z;

Klambda = exp((1j*t1)*lambda_s.')*exp((-1j*lambda_s)*t2.');
K = real(Klambda)/S;