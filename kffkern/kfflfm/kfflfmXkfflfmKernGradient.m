function [g1 ,g2] = kfflfmXkfflfmKernGradient(lfmKern1, lfmKern2, t1, t2, covGrad)

% KFFLFMXKFFLFMKERNGRADIENT Compute a cross gradient between two KFFLFM kernels.
%
%	Description:
%
%	[G1, G2] = KFFLFMXKFFLFMKERNGRADIENT(KFFSIMKERN1, KFFSIMKERN2, T, COVGRAD)
%	computes cross gradient of parameters of a cross kernel between two
%	sim kernels for the multiple output kernel.
%	 Returns:
%	  G1 - gradient of the parameters of the first kernel, for ordering
%	   see simKernExtractParam.
%	  G2 - gradient of the parameters of the second kernel, for ordering
%	   see simKernExtractParam.
%	 Arguments:
%	  KFFSIMKERN1 - the kernel structure associated with the first SIM
%	   kernel.
%	  KFFSIMKERN2 - the kernel structure associated with the second SIM
%	   kernel.
%	  T - inputs for which kernel is to be computed.
%	  COVGRAD - gradient of the objective function with respect to the
%	   elements of the cross kernel matrix.
%
%	[G1, G2] = KFFSIMXKFFSIMKERNGRADIENT(SIMKERN1, SIMKERN2, T1, T2, COVGRAD)
%	computes cross kernel terms between two SIM kernels for the multiple
%	output kernel.
%	 Returns:
%	  G1 - gradient of the parameters of the first kernel, for ordering
%	   see simKernExtractParam.
%	  G2 - gradient of the parameters of the second kernel, for ordering
%	   see simKernExtractParam.
%	 Arguments:
%	  SIMKERN1 - the kernel structure associated with the first SIM
%	   kernel.
%	  SIMKERN2 - the kernel structure associated with the second SIM
%	   kernel.
%	  T1 - row inputs for which kernel is to be computed.
%	  T2 - column inputs for which kernel is to be computed.
%	  COVGRAD - gradient of the objective function with respect to the
%	   elements of the cross kernel matrix.
%	
%	
%
%	See also
%	MULTIKERNPARAMINIT, MULTIKERNCOMPUTE, KFFLFMKERNPARAMINIT, KFFLFMKERNEXTRACTPARAM

%   Based on codes from N. Lawrence.
%	Copyright (c) 2018, Mauricio A. Alvarez, Cristian Guarnizo

if nargin < 4,
    covGrad = t2;
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
alpha(1) = lfmKern1.damper./(2*lfmKern1.mass);
alpha(2) = lfmKern2.damper./(2*lfmKern2.mass);
omega(1) = sqrt(lfmKern1.spring./lfmKern1.mass - alpha(1)*alpha(1));
omega(2) = sqrt(lfmKern2.spring./lfmKern2.mass - alpha(2)*alpha(2));

% using complex numbers
S = lfmKern1.S; 
Z = lfmKern1.Z;

sqInvW = sqrt(lfmKern1.inverseWidth);
lambda_s = sqInvW*Z'; %1xS
dlam_diw = 0.5/sqInvW*Z; %Sx1

%% wrt t1
elamt = exp(1j*t1*lambda_s); %NxS
temp = exp(-alpha(1)*t1);
esin = temp.*sin(omega(1)*t1); %Nx1
ecos = temp.*cos(omega(1)*t1); %Nx1
wm = omega(1)*lfmKern1.mass;

num = bsxfun(@minus,omega(1)*elamt - esin*(alpha(1)+1j*lambda_s),omega(1)*ecos);
den = wm*(alpha(1)^2+omega(1)^2-lambda_s.^2+2j*alpha(1)*lambda_s); %1xS
vd = bsxfun(@rdivide,num,den);

dvd_dm = vd*(-1/lfmKern1.mass);

temp = bsxfun(@times,1j*omega(1)*t1,elamt);
temp = bsxfun(@minus,temp,1j*esin);
temp = temp - bsxfun(@times,vd,wm*(-2*lambda_s+2j*alpha(1)));
dvd_dlam = bsxfun(@rdivide, temp, den);

temp = bsxfun(@plus,elamt-(t1.*ecos)*(alpha(1)+1j*lambda_s),omega(1)*(t1.*esin)-ecos);
temp = temp - bsxfun(@times,vd,lfmKern1.mass*(alpha(1)^2+3*omega(1)^2-lambda_s.^2+2j*alpha(1)*lambda_s));
dvd_dw = bsxfun(@rdivide, temp, den);

temp = bsxfun(@plus,(t1.*esin)*(alpha(1)+1j*lambda_s),omega(1)*(t1.*ecos)-esin);
temp = temp - bsxfun(@times,vd,2*wm*(alpha(1)+1j*lambda_s));
dvd_da = bsxfun(@rdivide, temp, den);

%% wrt t2
elamt = exp(-1j*t2*lambda_s); %NxS
temp = exp(-alpha(2)*t2);
esin = temp.*sin(omega(2)*t2); %Nx1
ecos = temp.*cos(omega(2)*t2); %Nx1
wm = omega(2)*lfmKern2.mass;

num = bsxfun(@minus,omega(2)*elamt - esin*(alpha(2)-1j*lambda_s),omega(2)*ecos); %NxS
den = wm*(alpha(2)^2+omega(2)^2-lambda_s.^2-2j*alpha(2)*lambda_s); %1xS
vdp = bsxfun(@rdivide,num,den);

dvdp_dm = vdp*(-1/lfmKern2.mass);

temp = bsxfun(@times,-1j*omega(2)*t2,elamt);
temp = bsxfun(@plus,temp,1j*esin);
temp = temp - bsxfun(@times,vdp,wm*(-2*lambda_s-2j*alpha(2)));
dvdp_dlam =  bsxfun(@rdivide, temp, den);

temp = bsxfun(@plus,elamt-(t2.*ecos)*(alpha(2)-1j*lambda_s),omega(2)*(t2.*esin)-ecos);
temp = temp - bsxfun(@times,vdp,lfmKern2.mass*(alpha(2)^2+3*omega(2)^2-lambda_s.^2-2j*alpha(2)*lambda_s));
dvdp_dw = bsxfun(@rdivide, temp, den);

temp = bsxfun(@plus,(t2.*esin)*(alpha(2)-1j*lambda_s),omega(2)*(t2.*ecos)-esin);
temp = temp - bsxfun(@times,vdp,2*wm*(alpha(2)-1j*lambda_s));
dvdp_da = bsxfun(@rdivide, temp, den);

%% Partial gradients
K = real(vd*vdp.');
dK_dm1 = real(dvd_dm*vdp.');
dK_dm2 = real(vd*dvdp_dm.');
dk_dinvWidth = real(dvd_dlam*bsxfun(@times,vdp.',dlam_diw)...
    +vd*bsxfun(@times,dvdp_dlam.',dlam_diw));
dK_dw1 = dvd_dw*vdp.';
dK_dw2 = vd*dvdp_dw.';
dK_da1 = real(dvd_da*vdp.');
dK_da2 = real(vd*dvdp_da.');

var2 = (lfmKern1.sensitivity*lfmKern2.sensitivity)/S;
%% Gradients
temp = sum(sum(covGrad.*dK_da1));
co_dw = sum(sum(covGrad.*dK_dw1));
dk_dm1 = var2*(sum(sum(covGrad.*dK_dm1)) + temp*(-alpha(1)/lfmKern1.mass)...
    + real(co_dw*(.5/omega(1)*(lfmKern1.damper^2/(2*lfmKern1.mass^3)-lfmKern1.spring/lfmKern1.mass^2))));
dk_ddamp1 = .5*var2/lfmKern1.mass*(temp + real(co_dw*(-alpha(1)/omega(1))));
dk_dspring1 = real(co_dw*var2*(.5/(lfmKern1.mass*omega(1))));

temp = sum(sum(covGrad.*dK_da2));
co_dw = sum(sum(covGrad.*dK_dw2));
dk_dm2 = var2*(sum(sum(covGrad.*dK_dm2)) + temp*(-alpha(2)/lfmKern2.mass)...
    + real(co_dw*(.5/omega(2)*(lfmKern2.damper^2/(2*lfmKern2.mass^3)-lfmKern2.spring/lfmKern2.mass^2))));
dk_ddamp2 = .5*var2/lfmKern2.mass*(temp + real(co_dw*(-alpha(2)/omega(2))));
dk_dspring2 = real(co_dw*var2*(.5/(lfmKern2.mass*omega(2))));

dk_dinvWidth = sum(sum(covGrad.*dk_dinvWidth))*var2;
temp = sum(sum(covGrad.*K));
dk_ds1 = lfmKern2.sensitivity/S * temp;
dk_ds2 = lfmKern1.sensitivity/S * temp;

% only pass the gradient with respect to the inverse width to one
% of the gradient vectors ... otherwise it is counted twice.
g1 = real([dk_dm1 dk_dspring1 dk_ddamp1 dk_dinvWidth dk_ds1]);
g2 = real([dk_dm2 dk_dspring2 dk_ddamp2 0 dk_ds2]);