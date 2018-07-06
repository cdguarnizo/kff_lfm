function gParam = dtcmgpLowerBoundGradients(model)

% IBPMULTIGPLOWERBOUNDGRADIENTS
% IBPMULTIGP
% MODIFICATIONS: Cristian Guarnizo, Mauricio Alvarez, 2014.

%[dLdKyy, dLdKuy, dLdKuu, dLdmu, gBeta] = spmultimodelLocalCovGradient2(model);
[dLdKyy, dLdKyu, dLdKuu, gBeta] = dtcmgpLowerBoundGradCovMat(model);
fhandle = str2func([model.betaTransform 'Transform']);
gBeta = gBeta.*fhandle(model.beta, 'gradfact');

if isfield(model, 'gamma') && ~isempty(model.gamma),
    gGamma = zeros(1, model.nlf);
    for i =1:model.nlf,
        gGamma(i) = trace(dLdKuu{i});
    end
    fhandle = str2func([model.gammaTransform 'Transform']);
    gGamma = gGamma.*fhandle(model.gamma, 'gradfact');
else
    gGamma = [];
end

gKernParam = kernGradient(model.kern, model.outX, model.latX, dLdKyy, dLdKyu, dLdKuu);

if isfield(model, 'fixInducing') && ~model.fixInducing,
    switch model.kern.type
        case 'ggglobal'
            gX = ggglobalKernGradient_X(model.kern, model.outX, model.latX, dLdKyu, dLdKuu);
        case 'lfmglobal'
            gX = lfmglobalKernGradient_X(model.kern, model.outX, model.latX, dLdKyu, dLdKuu);
        case 'simglobal'
            gX = simglobalKernGradient_X(model.kern, model.outX, model.latX, dLdKyu, dLdKuu);
    end
else
    gX = [];
end

gParam = [gKernParam gGamma gBeta gX];

if isfield(model, 'fix'),
    for i = 1:length(model.fix)
        gParam(model.fix(i).index) = 0;
    end
end