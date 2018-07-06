function [g1, g2] = kffsimXkffsimKernGradient(simKern1, simKern2, t1, t2, covGrad)

% KFFSIMXKFFSIMKERNGRADIENT Compute a cross gradient between two KFFSIM kernels.
%
%	Description:
%
%	[G1, G2] = KFFSIMXKFFSIMKERNGRADIENT(KFFSIMKERN1, KFFSIMKERN2, T, COVGRAD)
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
%	MULTIKERNPARAMINIT, MULTIKERNCOMPUTE, SIMKERNPARAMINIT, SIMKERNEXTRACTPARAM

%   Based on codes from N. Lawrence.
%	Copyright (c) 2018, Mauricio A. Alvarez, Cristian Guarnizo

arg{1}=t1;
if nargin < 5
  covGrad = t2;
  t2 = t1;
else
  arg{2}=t2;
end
if size(t1, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end
if simKern1.inverseWidth ~= simKern2.inverseWidth
  error('Kernels cannot be cross combined if they have different inverse widths.')
end

decay1 = simKern1.decay;
decay2 = simKern2.decay;
S = simKern1.S;
Z = simKern1.Z;
sqInvW = sqrt(simKern1.inverseWidth);
lambda_s = sqInvW*Z;
dlam_diw = 0.5/sqInvW*Z;

if isfield(simKern1, 'isNegativeS') && (simKern1.isNegativeS == true)
  C1 = simKern1.sensitivity;
  C2 = simKern2.sensitivity;
else
  C1 = sqrt(simKern1.variance);
  C2 = sqrt(simKern2.variance);
end

B = 1./(decay1 + 1i*lambda_s);
A = -B;
repB = repmat(B.', length(t1), 1);
ejlamt = exp(1i*t1*(lambda_s'));
eD1t1 = exp(-decay1*t1);
vd = repB.*ejlamt + eD1t1*(A.');

dvd_dlam = bsxfun(@times,(bsxfun(@times,t1,ejlamt) - vd),1j*B.');
dvd_dD1 = bsxfun(@times,bsxfun(@minus,t1.*eD1t1,vd),B.');

B = 1./(decay2 - 1i*lambda_s);
A = -B;
repB2 = repmat(B.', length(t2), 1);
ejlamt2 = exp(-1i*t2*(lambda_s'));
eD2t2= exp(-decay2*t2);
vdp = repB2.*ejlamt2 + eD2t2*(A.');

dvdp_dlam = bsxfun(@times,(bsxfun(@times,t2,ejlamt2) - vdp),-1j*B.');
dvdp_dD2 = bsxfun(@times,bsxfun(@minus,t2.*eD2t2,vdp),B.');

dK_dsigma = real(dvd_dlam*bsxfun(@times,vdp.',dlam_diw)...
    +vd*bsxfun(@times,dvdp_dlam.',dlam_diw));
dK_dD1 = real(dvd_dD1*vdp.');
dK_dD2 = real(vd*dvdp_dD2.');
K = real(vd*vdp.');

var2 = (C1*C2)/S;
dk_dD1 = sum(sum(covGrad.*dK_dD1))*var2;
dk_dD2 = sum(sum(covGrad.*dK_dD2))*var2;
dk_dinvWidth = sum(sum(covGrad.*dK_dsigma))*var2;
dk_dC1 = C2/S * sum(sum(covGrad.*K));
dk_dC2 = C1/S * sum(sum(covGrad.*K));

if isfield(simKern1, 'isNegativeS') && simKern1.isNegativeS
  dk_dSim1Variance = dk_dC1;
else
  dk_dSim1Variance = dk_dC1*0.5/C1;
end

if isfield(simKern2, 'isNegativeS') && simKern2.isNegativeS
  dk_dSim2Variance = dk_dC2;
else
  dk_dSim2Variance = dk_dC2*0.5/C2;
end

% only pass the gradient with respect to the inverse width to one
% of the gradient vectors ... otherwise it is counted twice.
g1 = real([dk_dD1 dk_dinvWidth dk_dSim1Variance]);
g2 = real([dk_dD2 0 dk_dSim2Variance]);
