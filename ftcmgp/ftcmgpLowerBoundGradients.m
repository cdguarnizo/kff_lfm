function gParam = ftcmgpLowerBoundGradients(model)

% FTCMGPLOWERBOUNDGRADIENTS

% Cristian Guarnizo, 2015

dLdKyy = model.alpha*model.alpha' - model.Kyyinv;
gBeta = [];
if model.includeNoise && isfield(model, 'beta'),
    gBeta = zeros(1,model.nout);
    temp = diag(dLdKyy);
    for d = 1:model.nout,
        gBeta(d) = (-1/model.beta(d)^2)*sum(temp(model.indX == d));
    end
    fhandle = str2func([model.betaTransform 'Transform']);
    gBeta = .5*gBeta.*fhandle(model.beta, 'gradfact');
end

gKernParam = kernGradient(model.kern, model.outX, mat2cell(dLdKyy,model.sizeX,model.sizeX));
gParam = [gKernParam gBeta];

if isfield(model, 'fix'),
    for i = 1:length(model.fix),
        gParam(model.fix(i).index) = 0;
    end
end