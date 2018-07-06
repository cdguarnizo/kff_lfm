function model = ftcmgpCreate( X, y, options )

% FTCMGPCREATE
% Inputs:
% X: Input data (time stamps and output indexation)
% y: Output data
% Options: Definition of flag variables and type of inference.
% FTCMGP

model.type = 'ftcmgp';

model.y = y;
model.nout = max(size(y));
model.ndim = size(X{1},2);
model.outX = X;
model.nlf = options.nlf;
model.approx = options.approx;
model.optimiser = options.optimiser;
model.kernType = options.kernType;

if isfield(options, 'scale') && ~isempty(options.scale)
    model.scale = options.scale;
else
    model.scale = ones(1, model.nout);
end
if isfield(options, 'bias') && ~isempty(options.bias)
    model.bias = options.bias;
else
    model.bias = zeros(1, model.nout);
end

model.tieIndices = options.tieOptions.tieIndices;
model.includeNoise = options.includeNoise;

% Set the dimension of each output in model, which will be useful for
% computing sizes and things like that

model.indX = [];
for k=1:model.nout,
    model.sizeX(k) = size(X{k},1);
    model.indX = [model.indX, k*ones(1,model.sizeX(k))];
end
model.N = sum(model.sizeX);

kern.type = [options.kernType 'global'];
kern.inputDimension = model.ndim;
kern.options = options.kern;
kern.options.nlf = options.nlf;
kern.options.nout = model.nout;
kern.options.approx = model.approx;

kern = kernParamInit(kern);
kern.template.output = kernCreate(X{1}, options.kernType);
if isfield(kern.template.output, 'isNegativeS'),
    kern.template.output.isNegativeS = 1;
end
kern.funcNames.computeOut = str2func([options.kernType 'KernCompute']);
kern.funcNames.computeCrossOut = str2func([options.kernType 'X' options.kernType 'KernCompute']);
kern.funcNames.gradientOut = str2func([options.kernType 'KernGradient']);
kern.funcNames.gradientCrossOut = str2func([options.kernType 'X' options.kernType 'KernGradient']);
kern.funcNames.extractOut = str2func([options.kernType 'KernExtractParam']);
model.kern = kern;
model.kern.isVarS = options.isVarS;
model.kernType = kern.type;
model.kern.paramGroups = speye(model.kern.nParams);
numParams = model.kern.nParams;

% Count number of parameters
model.nParams = numParams;

% Set up a mean function if one is given.
if isfield(options, 'meanFunction') && ~isempty(options.meanFunction)
    if isstruct(options.meanFunction)
        model.meanFunction = options.meanFunction;
    else
        if ~isempty(options.meanFunction)
            model.meanFunction = meanCreate(model.ndim, model.nout, X, y, options.meanFunctionOptions);
        end
    end
    model.nParams = model.nParams + model.meanFunction.nParams;
end

% Create noise models
switch model.approx,
    case 'ftc'
        if isfield(options, 'beta') && ~isempty(options.beta)
            if size(options.beta,2) == model.nout
                model.beta = options.beta;
            else
                model.beta = options.beta*ones(1,model.nout);
            end
            model.betaTransform =  optimiDefaultConstraint('positive');
            model.nParams = model.nParams + model.nout;
        end
    case {'fitc','pitc', 'dtcvar'}
        % In this structure is easier to put noise in the latent functions
        % at the top level
        if isfield(options, 'gamma') && ~isempty(options.gamma)
            if size(options.gamma,2) == model.nlf
                model.gamma = options.gamma;
            else
                model.gamma = options.gamma*ones(1,model.nlf);
            end
            model.gammaTransform =  optimiDefaultConstraint('positive');
            model.nParams = model.nParams + model.nlf;
        end
        if isfield(options, 'beta') && ~isempty(options.beta)
            if size(options.beta,2) == model.nout
                model.beta = options.beta;
            else
                model.beta = options.beta*ones(1,model.nout);
            end
            model.betaTransform =  optimiDefaultConstraint('positive');
            model.nParams = model.nParams + model.nout;
        end
end
params = modelExtractParam(model);
model = modelExpandParam(model, params);