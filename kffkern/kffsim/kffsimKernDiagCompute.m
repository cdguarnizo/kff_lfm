function [K, sK] = kffsimKernDiagCompute(simKern, t1)

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

decay1 = simKern.decay;
S = simKern.S;
Z = simKern.Z;
lambda_s = sqrt(simKern.inverseWidth)*Z;
B = 1./(decay1 + 1i*lambda_s);
A = -B;
repB = repmat(B.', length(t1), 1);
vd = repB.*exp(1i*t1*(lambda_s')) + exp(-decay1*t1)*(A.');

Klambda = sum(real(vd).^2+imag(vd).^2,2);
sK = 1/S*real(Klambda);

if isfield(simKern, 'isNegativeS') && (simKern.isNegativeS == true)
    K = simKern.sensitivity^2 * sK;
else
    K = sqrt(simKern.variance) *sK;
end