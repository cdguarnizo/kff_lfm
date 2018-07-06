function [f,g] = dtcmgpObjectiveGradient(params, model)

try
    model = dtcmgpExpandParam(model, params);
    model = dtcmgpComputeKernels(model);
    f = -dtcmgpLowerBound(model);
    g = -dtcmgpLowerBoundGradients(model);
catch
    f = inf;
    g = zeros(size(params));
end