function g = kffsimKernDiagGradient(simKern1, t1, covGrad)

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

decay1 = simKern1.decay;
S = simKern1.S;
Z = simKern1.Z;
sqInvW = sqrt(simKern1.inverseWidth);
lambda_s = sqInvW*Z;
dlam_diw = 0.5/sqInvW*Z;

if isfield(simKern1, 'isNegativeS') && (simKern1.isNegativeS == true)
  C1 = simKern1.sensitivity;
else
  C1 = sqrt(simKern1.variance);
end

B = 1./(decay1 + 1i*lambda_s);
A = -B;
repB = repmat(B.', length(t1), 1);
ejlamt = exp(1i*t1*(lambda_s'));
eD1t1 = exp(-decay1*t1);
vd = repB.*ejlamt + eD1t1*(A.');

dvd_dlam = bsxfun(@times,(bsxfun(@times,t1,ejlamt) - vd),1j*B.');
dvd_dD1 = bsxfun(@times,bsxfun(@minus,t1.*eD1t1,vd),B.');

vdp = conj(vd);

dK_dsigma = real(sum(dvd_dlam.*bsxfun(@times,vdp,dlam_diw.'),2)...
    +sum(vd.*bsxfun(@times,conj(dvd_dlam),dlam_diw.'),2));
dK_dD1 = real(sum(dvd_dD1.*vdp,2));
K = real(sum(vd.*vdp,2));

var2 = (C1*C1)/S;
dk_dD1 = sum(covGrad.*dK_dD1)*var2;
dk_dinvWidth = sum(covGrad.*dK_dsigma)*var2;
dk_dC1 = C1/S * sum(covGrad.*K);

if isfield(simKern1, 'isNegativeS') && simKern1.isNegativeS
  dk_dSim1Variance = dk_dC1;
else
  dk_dSim1Variance = dk_dC1*0.5/C1;
end

% only pass the gradient with respect to the inverse width to one
% of the gradient vectors ... otherwise it is counted twice.
g = real([2*dk_dD1 dk_dinvWidth 2*dk_dSim1Variance]);

