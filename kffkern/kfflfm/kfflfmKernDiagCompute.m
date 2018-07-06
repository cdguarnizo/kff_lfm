function k = kfflfmKernDiagCompute(Kern, t1)

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

% using complex numbers
S = Kern.S; 
Z = Kern.Z;

lambda_s = sqrt(Kern.inverseWidth)*Z';

elamt = exp(1j*t1*lambda_s); %NxS
temp = exp(-alpha*t1);
esin = temp.*sin(omega*t1); %Nx1
ecos = temp.*cos(omega*t1); %Nx1
wm = omega*Kern.mass;

num = bsxfun(@minus,omega*elamt - esin*(alpha+1j*lambda_s),omega*ecos);
den = wm*(alpha^2+omega^2-lambda_s.^2+2j*alpha*lambda_s); %1xS
vd = bsxfun(@rdivide,num,den);

num = bsxfun(@minus,omega*conj(elamt) - esin*(alpha-1j*lambda_s),omega*ecos);
den = wm*(alpha^2+omega^2-lambda_s.^2-2j*alpha*lambda_s); %1xS
vdp = bsxfun(@rdivide,num,den);

K0 = (Kern.sensitivity^2)/S;
k = K0*sum(vd.*vdp,2);
