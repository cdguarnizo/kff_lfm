function K = kfflfmXkffKernCompute(lfmKern1, lfmKern2, t1, t2)

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

if nargin < 4
    t2 = t1;
end
if size(t1, 2) > 1 || size(t2, 2) > 1
    error('Input can only have one column');
end
if lfmKern1.inverseWidth ~= lfmKern2.inverseWidth
    error('Kernels cannot be cross combined if they have different inverse widths.')
end

if lfmKern1.S ~= lfmKern2.S
    error('Kernels cannot be cross combined if they have different number of frequencies.')
end
    
% Parameters of the kernel
alpha = lfmKern1.damper./(2*lfmKern1.mass);
omega = sqrt(lfmKern1.spring./lfmKern1.mass - alpha*alpha);

% using complex numbers
S = lfmKern1.S; 
Z = lfmKern1.Z;
%Z = randn(S, 1);
lambda_s = sqrt(lfmKern1.inverseWidth)*Z;
C = computeC(alpha, omega, lambda_s);
repC = repmat(C.', length(t1), 1);
A = -C;
D = computeD(alpha, omega, lambda_s);
vd = repC.*exp(1i*t1*(lambda_s')) + (exp(-alpha*t1).*cos(omega*t1))*(A.') ...
    + (exp(-alpha*t1).*sin(omega*t1))*(D.');

Klambda = vd*exp((-1j*lambda_s)*t2.');
K0 = (lfmKern1.sensitivity)/(S*lfmKern1.mass*omega);
K = real(K0*Klambda);