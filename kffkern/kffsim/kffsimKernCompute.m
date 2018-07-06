function [k, sk, phi] = kffsimKernCompute(kern, t, t2)

% KFFSIMKERNCOMPUTE Compute the KFFSIM kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for the kff single input motif
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG t1 : the input matrix associated with the rows of the kernel.
% ARG t2 : the input matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
% RETURN sk : unscaled kernel matrix (i.e. only 0.5 times the sum of h's
% part).
%
% FORMAT
% DESC computes the kernel matrix for the kff single input motif
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG t : input data matrix in the form of a design matrix.
% RETURN k : the kernel matrix computed at the given points.
% RETURN sk : unscaled kernel matrix (i.e. only 0.5 times the sum of h's
% part).
%
% SEEALSO : kffsimKernParamInit, kernCompute, kernCreate, simKernDiagCompute
%
% COPYRIGHT :  Mauricio Alvarez, 2017


% KERN

flag = 0;
if nargin < 3
  t2 = t;
  flag= 1;
end
if size(t, 2) > 1 || size(t2, 2) > 1
  error('Input can only have one column');
end

decay = kern.decay;
lengthscale2 = 2/kern.inverseWidth;
S = kern.S; 
Z = kern.Z;
lambda_s = sqrt(2/lengthscale2)*Z;

B = 1./(decay + 1i*lambda_s);
A = -B;
repB = repmat(B.', length(t), 1);
vd = repB.*exp(1i*t*(lambda_s')) + exp(-decay*t)*(A.');

if flag,
    Klambda = vd*(vd');
else
    B = 1./(decay - 1i*lambda_s);
    A = -B;
    repB2 = repmat(B.', length(t2), 1);
    vdp = repB2.*exp(-1i*t2*(lambda_s')) + exp(-decay*t2)*(A.');
    
    Klambda = vd*vdp.';
end
    
sk = Klambda/S;
if isfield(kern, 'isNegativeS') && kern.isNegativeS
  K0 = kern.sensitivity^2;  
else
  K0 = kern.variance;
end
k = real(K0*sk);
if nargout > 2,
    phi = sqrt(K0/S)*vd;
end
