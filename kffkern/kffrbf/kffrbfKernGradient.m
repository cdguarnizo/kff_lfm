function g = kffrbfKernGradient(Kern, t, covGrad)

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

% using complex numbers
S = Kern.S;
Z = Kern.Z;
sqInvW = sqrt(Kern.inverseWidth);
lambda_s = sqInvW*Z; %1xS
dlam_diw = 0.5/sqInvW*Z; %Sx1

kd = exp((1j*t)*lambda_s.');
kdp = conj(kd);

dkd_dlam =  bsxfun(@times,1j*t,kd);
dkdp_dlam =  conj(dkd_dlam);

dk_dinvWidth = real(dkd_dlam*bsxfun(@times,kdp.',dlam_diw)...
    +kd*bsxfun(@times,dkdp_dlam.',dlam_diw));

g(1) = sum(sum(covGrad.*dk_dinvWidth))/S;
g(2) = 0;