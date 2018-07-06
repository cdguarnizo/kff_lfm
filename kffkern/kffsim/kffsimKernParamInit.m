function kern = kffsimKernParamInit(kern)
% KFFSIMKERNPARAMINIT Initializes the kernel Fourier feature sim. It
% consists of the SIM kernel plus the parameter S of the number of samples
% in the Gaussian distribution
%
% COPYRIGHT : Mauricio A. Alvarez, 2017

% KERN

kern = simKernParamInit(kern);
kern.type = 'kffsim';
if isfield(kern, 'options') && isfield(kern.options, 'S') 
    kern.S = kern.options.S;    
else
    kern.S = 10;
end
if isfield(kern, 'options') && isfield(kern.options, 'Z') 
    kern.Z = kern.options.Z;    
else
    kern.Z = randn(kern.S, 1); % Create the basis functions in the kernel Fourier feature way
end