function kern = kfflfmKernParamInit(kern)
% KFFLFMKERNPARAMINIT Initializes the kernel Fourier feature lfm. It
% consists of the LFM kernel plus the parameter S of the number of samples
% in the Gaussian distribution
%
% COPYRIGHT : Mauricio A. Alvarez, 2017

% KERN

kern = lfmKernParamInit(kern);
kern.type = 'kfflfm';
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
