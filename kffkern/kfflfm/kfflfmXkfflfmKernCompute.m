function K = kfflfmXkfflfmKernCompute(lfmKern1, lfmKern2, t1, t2)

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
if lfmKern1.inverseWidth ~= lfmKern2.inverseWidth
    error('Kernels cannot be cross combined if they have different inverse widths.')
end

if lfmKern1.S ~= lfmKern2.S
    error('Kernels cannot be cross combined if they have different number of frequencies.')
end
    
% Get length scale out.
lengthscale2 = 2/lfmKern1.inverseWidth; % This is what we call \ell_q2 for the SE kernel in LFM

% Parameters of the kernel
alpha(1) = lfmKern1.damper./(2*lfmKern1.mass);
alpha(2) = lfmKern2.damper./(2*lfmKern2.mass);
omega(1) = sqrt(lfmKern1.spring./lfmKern1.mass - alpha(1)*alpha(1));
omega(2) = sqrt(lfmKern2.spring./lfmKern2.mass - alpha(2)*alpha(2));

% Sample the frequencies for the Monte Carlo integration
% S = lfmKern1.S; 
% lambda_s = sqrt(2/lengthscale2)*randn(S, 1);
% Klambda = zeros(length(t1), length(t2));
% for s=1:S
%     zt  = computez(t1, alpha(1), omega(1), lambda_s(s));
%     ztp = computez(t2, alpha(2), omega(2), lambda_s(s));
%     Klambda = Klambda + zt*(ztp.');    
% end
% Using matrix computations
% S = lfmKern1.S; 
% lambda_s = sqrt(2/lengthscale2)*randn(S, 1);
% tt = repmat(t1, 1, S);
% ttp = repmat(t2, 1, S);
% Lambda_s_tt = repmat(lambda_s', length(t1), 1);
% Lambda_s_ttp = repmat(lambda_s', length(t2), 1);
% [Z1t, Z2t] = computeZmatrix(tt, alpha(1), omega(1), Lambda_s_tt);
% [Z1tp, Z2tp] = computeZmatrix(ttp, alpha(2), omega(2), Lambda_s_ttp);
% Klambda = Z1t*(Z1tp.') + Z2t*(Z2tp.');
% Klambda = Klambda/S;
% sK = (exp(-alpha(1)*t1)*(exp(-alpha(2)*t2)).').*Klambda;
% K0 = (lfmKern1.sensitivity*lfmKern2.sensitivity)/(lfmKern1.mass*lfmKern2.mass*prod(omega));
% K = K0*sK;
% using complex numbers
S = lfmKern1.S; 
Z = lfmKern1.Z;
%Z = randn(S, 1);
lambda_s = sqrt(2/lengthscale2)*Z;
C = computeC(alpha(1), omega(1), lambda_s);
repC = repmat(C.', length(t1), 1);
A = -C;
D = computeD(alpha(1), omega(1), lambda_s);
vd = repC.*exp(1i*t1*(lambda_s')) + (exp(-alpha(1)*t1).*cos(omega(1)*t1))*(A.') ...
    + (exp(-alpha(1)*t1).*sin(omega(1)*t1))*(D.');
C2 = computeC(alpha(2), omega(2), -lambda_s);
repC2 = repmat(C2.', length(t2), 1);
A2 = -C2;
D2 = computeD(alpha(2), omega(2), -lambda_s);
vdp = repC2.*exp(-1i*t2*(lambda_s')) + (exp(-alpha(2)*t2).*cos(omega(2)*t2))*(A2.') ...
    + (exp(-alpha(2)*t2).*sin(omega(2)*t2))*(D2.');
Klambda = vd*(vdp.');
sK = Klambda/S;
K0 = (lfmKern1.sensitivity*lfmKern2.sensitivity)/(lfmKern1.mass*lfmKern2.mass*prod(omega));
K = real(K0*sK);


% if isreal(omega)
%     % Precomputations to increase speed
%     gamma1 = alpha(1) + j*omega(1);
%     gamma2 = alpha(2) + j*omega(2);
%     preGamma(1) = gamma1 + gamma2;
%     preGamma(2) = conj(gamma1) + gamma2;
%     preConst = 1./preGamma;
%     preExp1 = exp(-gamma1*t1);
%     preExp2 = exp(-gamma2*t2);
%     % Actual computation of the kernel
%     sK = real(lfmComputeH3(gamma1, gamma2, sigma2, t1,t2,preConst, 0, 1) + ...
%         lfmComputeH3(gamma2, gamma1, sigma2, t2,t1,preConst(2) - preConst(1), 0, 0).' + ...
%         lfmComputeH4(gamma1, gamma2, sigma2, t1, preGamma, preExp2, 0, 1  ) + ...
%         lfmComputeH4(gamma2, gamma1, sigma2, t2, preGamma, preExp1,0, 0 ).');
%     if lfmKern1.isNormalised
%         K0 = (lfmKern1.sensitivity*lfmKern2.sensitivity)/(4*sqrt(2)*lfmKern1.mass*lfmKern2.mass*prod(omega));
%     else
%         K0 = (sigma*sqrt(pi)*lfmKern1.sensitivity*lfmKern2.sensitivity)/(4*lfmKern1.mass*lfmKern2.mass*prod(omega));
%     end
%     K = K0*sK;    
% else
%     % Precomputations to increase the speed
%     preExp1 = zeros(length(t1),2);
%     preExp2 = zeros(length(t2),2);
%     gamma1_p = alpha(1) + j*omega(1);
%     gamma1_m = alpha(1) - j*omega(1);
%     gamma2_p = alpha(2) + j*omega(2);
%     gamma2_m = alpha(2) - j*omega(2);
%     preGamma(1) = gamma1_p + gamma2_p;
%     preGamma(2) = gamma1_p + gamma2_m;
%     preGamma(3) = gamma1_m + gamma2_p;
%     preGamma(4) = gamma1_m + gamma2_m;
%     preConst = 1./preGamma;
%     preFactors(1) = preConst(2) - preConst(1);
%     preFactors(2) = preConst(3) - preConst(4);
%     preFactors(3) = preConst(3) - preConst(1);
%     preFactors(4) = preConst(2) - preConst(4);
%     preExp1(:,1) = exp(-gamma1_p*t1);
%     preExp1(:,2) = exp(-gamma1_m*t1);
%     preExp2(:,1) = exp(-gamma2_p*t2);
%     preExp2(:,2) = exp(-gamma2_m*t2);
%     % Actual computation of the kernel
%     sK = (  lfmComputeH3(gamma1_p, gamma1_m, sigma2, t1,t2,preFactors([1 2]), 1) + ...
%         lfmComputeH3(gamma2_p, gamma2_m, sigma2, t2,t1,preFactors([3 4]), 1).' + ...
%         lfmComputeH4(gamma1_p, gamma1_m, sigma2, t1, preGamma([1 2 4 3]), preExp2, 1 ) + ...
%         lfmComputeH4(gamma2_p, gamma2_m, sigma2, t2, preGamma([1 3 4 2]), preExp1, 1 ).');
%     if lfmKern1.isNormalised
%         K0 = (lfmKern1.sensitivity*lfmKern2.sensitivity)/(8*sqrt(2)*lfmKern1.mass*lfmKern2.mass*prod(omega));
%     else
%         K0 = (sigma*sqrt(pi)*lfmKern1.sensitivity*lfmKern2.sensitivity)/(8*lfmKern1.mass*lfmKern2.mass*prod(omega));
%     end
%     K = K0*sK;
% end
