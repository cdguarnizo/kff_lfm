function [K, sK] = kffsimXkffsimKernCompute(simKern1, simKern2, t1, t2)

% KFFSIMXKFFSIMKERNCOMPUTE Compute a cross kernel between two KFF SIM kernels.
% FORMAT
% DESC computes cross kernel terms between two sim kernels for
% the multiple output kernel.
% ARG simKern1 : the kernel structure associated with the first sim
% kernel.
% ARG simKern2 : the kernel structure associated with the second sim
% kernel.
% ARG t : inputs for which kernel is to be computed.
% RETURN K : block of values from kernel matrix.
%
% FORMAT
% DESC computes cross kernel terms between two KFF sim kernels for
% the multiple output kernel.
% ARG simKern1 : the kernel structure associated with the first sim
% kernel.
% ARG simKern2 : the kernel structure associated with the second sim
% kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% RETURN K : block of values from kernel matrix.
%
% SEEALSO : kffsimKernParamInit, kffsimKernCompute, kffsimKernParamInit
%
% COPYRIGHT : Mauricio Alvarez, 2017

% KERN

if nargin < 4
    t2 = t1;
end
if size(t1, 2) > 1 || size(t2, 2) > 1
    error('Input can only have one column');
end
if simKern1.inverseWidth ~= simKern2.inverseWidth
    error('Kernels cannot be cross combined if they have different inverse widths.')
end

decay1 = simKern1.decay;
decay2 = simKern2.decay;
S = simKern1.S;
Z = simKern1.Z;
lambda_s = sqrt(simKern1.inverseWidth)*Z;
B = 1./(decay1 + 1i*lambda_s);
A = -B;
repB = repmat(B.', length(t1), 1);
vd = repB.*exp(1i*t1*(lambda_s')) + exp(-decay1*t1)*(A.');
B = 1./(decay2 - 1i*lambda_s);
A = -B;
repB2 = repmat(B.', length(t2), 1);
vdp = repB2.*exp(-1i*t2*(lambda_s')) + exp(-decay2*t2)*(A.');
Klambda = vd*(vdp.');
sK = 1/S*real(Klambda);

if isfield(simKern1, 'isNegativeS') && (simKern1.isNegativeS == true)
    K = simKern1.sensitivity * sK;
else
    K = sqrt(simKern1.variance) *sK;
end
if isfield(simKern2, 'isNegativeS') && (simKern2.isNegativeS == true)
    K = simKern2.sensitivity * K;
else
    K = sqrt(simKern2.variance) * K;
end

