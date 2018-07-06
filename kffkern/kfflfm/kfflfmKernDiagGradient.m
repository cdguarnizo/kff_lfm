function g = kfflfmKernDiagGradient(Kern, t1, covGrad)

% KFFLFMKERNCOMPUTE Compute the KFFLFM kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for the KFF latent force model
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG t1 : the input matrix associated with the rows of the kernel.
% ARG t2 : the input matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel matrix for the KFF latent force model
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG t : input data matrix in the form of a design matrix.
% RETURN k : the kernel matrix computed at the given points.
%
% SEEALSO : lfmKernParamInit, kernCompute, kernCreate, lfmKernDiagCompute
%
% COPYRIGHT : Mauricio A. Alvarez, 2017

% KERN


% Parameters of the kernel
alpha = Kern.damper./(2*Kern.mass);
omega = sqrt(Kern.spring./Kern.mass - alpha*alpha);

S = Kern.S; 
Z = Kern.Z;

sqInvW = sqrt(Kern.inverseWidth);
lambda_s = sqInvW*Z'; %1xS
dlam_diw = 0.5/sqInvW*Z; %Sx1

%% wrt t1
elamt = exp(1j*t1*lambda_s); %NxS
temp = exp(-alpha*t1);
esin = temp.*sin(omega*t1); %Nx1
ecos = temp.*cos(omega*t1); %Nx1
wm = omega*Kern.mass;

num = bsxfun(@minus,omega*elamt - esin*(alpha+1j*lambda_s),omega*ecos);
den = wm*(alpha^2+omega^2-lambda_s.^2+2j*alpha*lambda_s); %1xS
vd = bsxfun(@rdivide,num,den);

dvd_dm = vd*(-1/Kern.mass);

temp = bsxfun(@times,1j*omega*t1,elamt);
temp = bsxfun(@minus,temp,1j*esin);
temp = temp - bsxfun(@times,vd,wm*(-2*lambda_s+2j*alpha));
dvd_dlam = bsxfun(@rdivide, temp, den);

temp = bsxfun(@plus,elamt-(t1.*ecos)*(alpha+1j*lambda_s),omega*(t1.*esin)-ecos);
temp = temp - bsxfun(@times,vd,Kern.mass*(alpha^2+3*omega^2-lambda_s.^2+2j*alpha*lambda_s));
dvd_dw = bsxfun(@rdivide, temp, den);

temp = bsxfun(@plus,(t1.*esin)*(alpha+1j*lambda_s),omega*(t1.*ecos)-esin);
temp = temp - bsxfun(@times,vd,2*wm*(alpha+1j*lambda_s));
dvd_da = bsxfun(@rdivide, temp, den);

%%
elamt = conj(elamt);
num = bsxfun(@minus,omega*elamt - esin*(alpha-1j*lambda_s),omega*ecos);
den = wm*(alpha^2+omega^2-lambda_s.^2-2j*alpha*lambda_s); %1xS
vdp = bsxfun(@rdivide,num,den);

% temp = bsxfun(@times,-1j*omega*t1,elamt);
% temp = bsxfun(@minus,temp,1j*esin);
% temp = temp - bsxfun(@times,vdp,wm*(-2*lambda_s-2j*alpha));
% dvdp_dlam = bsxfun(@rdivide, temp, den);

%% Partial gradients
K = sum(vd.*vdp,2);
dK_dm1 = real(sum(dvd_dm.*vdp,2));
dk_dinvWidth = real(sum(bsxfun(@times,dvd_dlam,dlam_diw.').*vdp,2));%...
%    +sum(vd.*bsxfun(@times,dvdp_dlam,dlam_diw.'),2));
dK_dw1 = sum(dvd_dw.*vdp,2);
dK_da1 = real(sum(dvd_da.*vdp,2));

var2 = (Kern.sensitivity^2)/S;
%% Gradients
temp = sum(covGrad.*dK_da1);
co_dw = sum(covGrad.*dK_dw1);
dk_dm1 = var2*(sum(covGrad.*dK_dm1) + temp*(-alpha/Kern.mass)...
    + real(co_dw*(.5/omega*(Kern.damper^2/(2*Kern.mass^3)-Kern.spring/Kern.mass^2))));
dk_ddamp1 = .5*var2/Kern.mass*(temp + real(co_dw*(-alpha/omega)));
dk_dspring1 = real(co_dw*var2*(.5/(Kern.mass*omega)));

dk_dinvWidth = sum(covGrad.*dk_dinvWidth)*var2;
temp = sum(covGrad.*K);
dk_ds1 = Kern.sensitivity/S * temp;

% only pass the gradient with respect to the inverse width to one
% of the gradient vectors ... otherwise it is counted twice.
g = 2.*[dk_dm1 dk_dspring1 dk_ddamp1 dk_dinvWidth dk_ds1];