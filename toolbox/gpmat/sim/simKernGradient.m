function g = simKernGradient(kern, t, varargin)

% SIMKERNGRADIENT Gradient of SIM kernel's parameters.
%
%	Description:
%
%	G = SIMKERNGRADIENT(KERN, T, PARTIAL) computes the gradient of
%	functions with respect to the single input motif kernel's
%	parameters. As well as the kernel structure and the input positions,
%	the user provides a matrix PARTIAL which gives the partial
%	derivatives of the function with respect to the relevant elements of
%	the kernel matrix.
%	 Returns:
%	  G - gradients of the function of interest with respect to the
%	   kernel parameters. The ordering of the vector should match that
%	   provided by the function kernExtractParam.
%	 Arguments:
%	  KERN - the kernel structure for which the gradients are being
%	   computed.
%	  T - the input locations for which the gradients are being
%	   computed.
%	  PARTIAL - matrix of partial derivatives of the function of
%	   interest with respect to the kernel matrix. The argument takes the
%	   form of a square matrix of dimension  numData, where numData is
%	   the number of rows in X.
%
%	G = SIMKERNGRADIENT(KERN, T1, T2, PARTIAL) computes the derivatives
%	as above, but input locations are now provided in two matrices
%	associated with rows and columns of the kernel matrix.
%	 Returns:
%	  G - gradients of the function of interest with respect to the
%	   kernel parameters.
%	 Arguments:
%	  KERN - the kernel structure for which the gradients are being
%	   computed.
%	  T1 - the input locations associated with the rows of the kernel
%	   matrix.
%	  T2 - the input locations associated with the columns of the kernel
%	   matrix.
%	  PARTIAL - matrix of partial derivatives of the function of
%	   interest with respect to the kernel matrix. The matrix should have
%	   the same number of rows as X1 and the same number of columns as X2
%	   has rows.
%	
%
%	See also
%	% SEEALSO SIMKERNPARAMINIT, KERNGRADIENT, SIMKERNDIAGGRADIENT, KERNGRADX


%	Copyright (c) 2006 Neil D. Lawrence


if length(varargin)<2
  t2 = t;
else
  t2 = varargin{1};
end

[g1, g2] = simXsimKernGradient(kern, kern, t, t2, varargin{end});

g = g1 + g2;

if isfield(kern, 'gaussianInitial') && kern.gaussianInitial,
  dim1 = size(t, 1);
  dim2 = size(t2, 1);
  t1Mat = t(:, ones(1, dim2));
  t2Mat = t2(:, ones(1, dim1))';

  % Gradient with respect to the decay
  g(1) = g(1) + sum(sum(-kern.initialVariance*(t1Mat + t2Mat).*exp(-kern.decay*(t1Mat + t2Mat)) .* varargin{end}));
  g(end+1) = sum(sum(exp(-kern.decay*(t1Mat + t2Mat)) .* varargin{end}));
end
