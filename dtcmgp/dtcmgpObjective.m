function f = dtcmgpObjective(params, model)

% DTCMGPOBJECTIVE Wrapper function for MODELOPTIMISE objective.

% DTCMGP

try
    model = dtcmgpExpandParam(model, params);
    model = dtcmgpComputeKernels(model); %Update Kff Kfu and Kuu
    f = -dtcmgpLowerBound(model);
catch
    f = inf;
end