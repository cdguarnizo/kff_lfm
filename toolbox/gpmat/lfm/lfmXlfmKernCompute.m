function K = lfmXlfmKernCompute(lfmKern1, lfmKern2, t1, t2)

% LFMXLFMKERNCOMPUTE Compute a cross kernel between two LFM kernels.
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
% DESC computes cross kernel terms between two LFM kernels for
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
% COPYRIGHT : David Luengo, 2007, 2008, Mauricio Alvarez, 2008
%
% MODIFICATIONS : Neil D. Lawrence, 2007, 2008,
%

% KERN

if nargin < 4
    t2 = t1;
end
if size(t1, 2) > 1 || size(t2, 2) > 1
    error('Input can only have one column');
end
if lfmKern1.inverseWidth ~= lfmKern2.inverseWidth
    error('Kernels cannot be cross combined if they have different inverse widths.')
end
    
% Get length scale out.
sigma2 = 2/lfmKern1.inverseWidth;
sigma = sqrt(sigma2);

% Parameters of the kernel
alpha(1) = lfmKern1.damper./(2*lfmKern1.mass);
alpha(2) = lfmKern2.damper./(2*lfmKern2.mass);
omega(1) = sqrt(lfmKern1.spring./lfmKern1.mass - alpha(1)*alpha(1));
omega(2) = sqrt(lfmKern2.spring./lfmKern2.mass - alpha(2)*alpha(2));

% Creation of the time matrices


if isreal(omega)
    % Precomputations to increase speed
    gamma1 = alpha(1) + j*omega(1);
    gamma2 = alpha(2) + j*omega(2);
    preGamma(1) = gamma1 + gamma2;
    preGamma(2) = conj(gamma1) + gamma2;
    preConst = 1./preGamma;
    preExp1 = exp(-gamma1*t1);
    preExp2 = exp(-gamma2*t2);
    % Actual computation of the kernel
    sK = real(lfmComputeH3(gamma1, gamma2, sigma2, t1,t2,preConst, 0, 1) + ...
        lfmComputeH3(gamma2, gamma1, sigma2, t2,t1,preConst(2) - preConst(1), 0, 0).' + ...
        lfmComputeH4(gamma1, gamma2, sigma2, t1, preGamma, preExp2, 0, 1  ) + ...
        lfmComputeH4(gamma2, gamma1, sigma2, t2, preGamma, preExp1,0, 0 ).');
    if lfmKern1.isNormalised
        K0 = (lfmKern1.sensitivity*lfmKern2.sensitivity)/(4*sqrt(2)*lfmKern1.mass*lfmKern2.mass*prod(omega));
    else
        K0 = (sigma*sqrt(pi)*lfmKern1.sensitivity*lfmKern2.sensitivity)/(4*lfmKern1.mass*lfmKern2.mass*prod(omega));
    end
    K = K0*sK;    
else
    % Precomputations to increase the speed
    preExp1 = zeros(length(t1),2);
    preExp2 = zeros(length(t2),2);
    gamma1_p = alpha(1) + j*omega(1);
    gamma1_m = alpha(1) - j*omega(1);
    gamma2_p = alpha(2) + j*omega(2);
    gamma2_m = alpha(2) - j*omega(2);
    preGamma(1) = gamma1_p + gamma2_p;
    preGamma(2) = gamma1_p + gamma2_m;
    preGamma(3) = gamma1_m + gamma2_p;
    preGamma(4) = gamma1_m + gamma2_m;
    preConst = 1./preGamma;
    preFactors(1) = preConst(2) - preConst(1);
    preFactors(2) = preConst(3) - preConst(4);
    preFactors(3) = preConst(3) - preConst(1);
    preFactors(4) = preConst(2) - preConst(4);
    preExp1(:,1) = exp(-gamma1_p*t1);
    preExp1(:,2) = exp(-gamma1_m*t1);
    preExp2(:,1) = exp(-gamma2_p*t2);
    preExp2(:,2) = exp(-gamma2_m*t2);
    % Actual computation of the kernel
    sK = (  lfmComputeH3(gamma1_p, gamma1_m, sigma2, t1,t2,preFactors([1 2]), 1) + ...
        lfmComputeH3(gamma2_p, gamma2_m, sigma2, t2,t1,preFactors([3 4]), 1).' + ...
        lfmComputeH4(gamma1_p, gamma1_m, sigma2, t1, preGamma([1 2 4 3]), preExp2, 1 ) + ...
        lfmComputeH4(gamma2_p, gamma2_m, sigma2, t2, preGamma([1 3 4 2]), preExp1, 1 ).');
    if lfmKern1.isNormalised
        K0 = (lfmKern1.sensitivity*lfmKern2.sensitivity)/(8*sqrt(2)*lfmKern1.mass*lfmKern2.mass*prod(omega));
    else
        K0 = (sigma*sqrt(pi)*lfmKern1.sensitivity*lfmKern2.sensitivity)/(8*lfmKern1.mass*lfmKern2.mass*prod(omega));
    end
    K = K0*sK;
end
