function g = dtcmgpGradient(params, model)

% DTCMGPGRADIENT 

% DTCMGP

try
    model = dtcmgpExpandParam(model, params);
    model = dtcmgpComputeKernels(model); %Update Kff Kfu and Kuu
    g = -dtcmgpLowerBoundGradients(model);
catch
    g = zeros(size(params));
end