function  [param, names] = dtcmgpExtractParam(model)

% IBPMULTIGPEXTRACTPARAM Extract the parameters of an IBPMULTIGP struct.

% IBPMULTIGP


if nargout > 1
    [paramKern, namesKern] = kernExtractParam(model.kern);    
else
    paramKern = kernExtractParam(model.kern);
end

if isfield(model, 'gamma') && ~isempty(model.gamma)
    fhandle = str2func([model.gammaTransform 'Transform']);
    gammaParams = fhandle(model.gamma, 'xtoa');    
    if nargout>1
        gammaParamNames = cell(model.nlf,1);
        for i = 1:length(gammaParams)
            gammaParamNames{i} = ['Gamma ' num2str(i)];
        end
    end
else
    gammaParamNames = {};
    gammaParams =[];
end

if isfield(model, 'beta') && ~isempty(model.beta)
    fhandle = str2func([model.betaTransform 'Transform']);
    betaParams = fhandle(model.beta, 'xtoa');    
    if nargout>1
        betaParamNames = cell(model.nout,1);
        for i = 1:length(betaParams)
            betaParamNames{i} = ['Beta ' num2str(i)];
        end
    end
else
    betaParamNames = {};
    betaParams =[];
end

if isfield(model, 'fixInducing') && ~model.fixInducing,
    %We assume that all latent functions have the same number of inducing
    %points, each point in R^p
    width = size(model.latX{1},1)*size(model.latX{1},2);
    XParams = zeros(1, model.nlf*width);
    XParamNames = cell(model.nlf*width, 1);
    for k=1:model.nlf,
        XParams(1 + (k-1)*width: k*width) = model.latX{k}(:)';
        for l=1:width,
            XParamNames{(k-1)*width + l} = ['X latent ' num2str(k)];
        end
    end
else
    XParamNames = {};
    XParams = [];
end

param = [paramKern gammaParams betaParams XParams];

% Fix the value of the parameters

if isfield(model, 'fix')
    for i = 1:length(model.fix)
        param(model.fix(i).index) = model.fix(i).value;
    end
end

if nargout > 1
    names = {namesKern{:}, gammaParamNames{:}, betaParamNames{:}, XParamNames{:}};
end