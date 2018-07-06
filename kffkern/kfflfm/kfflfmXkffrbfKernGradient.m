function [g1,g2] = kfflfmXkffrbfKernGradient(lfmKern, rbfKern, t1, t2, covGrad)

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
if lfmKern.inverseWidth ~= rbfKern.inverseWidth
    error('Kernels cannot be cross combined if they have different inverse widths.')
end

if lfmKern.S ~= rbfKern.S
    error('Kernels cannot be cross combined if they have different number of frequencies.')
end
    
% Parameters of the kernel
alpha = lfmKern.damper./(2*lfmKern.mass);
omega = sqrt(lfmKern.spring./lfmKern.mass - alpha*alpha);

% using complex numbers
S = lfmKern.S; 
Z = lfmKern.Z;
%Z = randn(S, 1);
sqInvW = sqrt(rbfKern.inverseWidth);
lambda_s = sqInvW*Z.';
dlam_diw = 0.5/sqInvW*Z; %Sx1

%% wrt t1
elamt = exp(1j*t1*lambda_s); %NxS
temp = exp(-alpha*t1);
esin = temp.*sin(omega*t1); %Nx1
ecos = temp.*cos(omega*t1); %Nx1
wm = omega*lfmKern.mass;

num = bsxfun(@minus,omega*elamt - esin*(alpha+1j*lambda_s),omega*ecos);
den = wm*(alpha^2+omega^2-lambda_s.^2+2j*alpha*lambda_s); %1xS
vd = bsxfun(@rdivide,num,den); %NxS

dvd_dm = vd*(-1/lfmKern.mass);

temp = bsxfun(@times,1j*omega*t1,elamt);
temp = bsxfun(@minus,temp,1j*esin);
temp = temp - bsxfun(@times,vd,wm*(-2*lambda_s+2j*alpha));
dvd_dlam = bsxfun(@rdivide, temp, den);

temp = bsxfun(@plus,elamt-(t1.*ecos)*(alpha+1j*lambda_s),omega*(t1.*esin)-ecos);
temp = temp - bsxfun(@times,vd,lfmKern.mass*(alpha^2+3*omega^2-lambda_s.^2+2j*alpha*lambda_s));
dvd_dw = bsxfun(@rdivide, temp, den);

temp = bsxfun(@plus,(t1.*esin)*(alpha+1j*lambda_s),omega*(t1.*ecos)-esin);
temp = temp - bsxfun(@times,vd,2*wm*(alpha+1j*lambda_s));
dvd_da = bsxfun(@rdivide, temp, den);

%% wrt t2
vdp = exp(t2*(-1j*lambda_s)); %NxS

dvdp_dlam = bsxfun(@times,-1j*t2,vdp); %NxS

%% Partial gradients
K = real(vd*vdp.');
dK_dm1 = real(dvd_dm*vdp.');
dk_dinvWidth = real(dvd_dlam*bsxfun(@times,vdp.',dlam_diw)...
    +vd*bsxfun(@times,dvdp_dlam.',dlam_diw));
dK_dw1 = dvd_dw*vdp.';
dK_da1 = real(dvd_da*vdp.');

var2 = lfmKern.sensitivity/S;
%% Gradients
temp = sum(sum(covGrad.*dK_da1));
co_dw = sum(sum(covGrad.*dK_dw1));
dk_dm1 = var2*(sum(sum(covGrad.*dK_dm1)) + temp*(-alpha/lfmKern.mass)...
    + real(co_dw*(.5/omega*(lfmKern.damper^2/(2*lfmKern.mass^3)-lfmKern.spring/lfmKern.mass^2))));
dk_ddamp1 = .5*var2/lfmKern.mass*(temp + real(co_dw*(-alpha/omega)));
dk_dspring1 = real(co_dw*var2*(.5/(lfmKern.mass*omega)));

dk_dinvWidth = sum(sum(covGrad.*dk_dinvWidth))*var2;
temp = sum(sum(covGrad.*K));
dk_ds1 = temp/S;

% only pass the gradient with respect to the inverse width to one
% of the gradient vectors ... otherwise it is counted twice.
g1 = real([dk_dm1 dk_dspring1 dk_ddamp1 dk_dinvWidth dk_ds1]);
g2 = real([dk_dinvWidth 0]);