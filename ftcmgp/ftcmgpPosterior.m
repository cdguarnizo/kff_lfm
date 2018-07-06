function [ymean, yval] = ftcmultigpPosterior(model, Xtest)

% FTCMULTIGPPOSTERIOR
% FTCMULTIGP

ymean = cell(model.nout, 1);
yvar = cell(model.nout, 1);

fhandle = str2func([model.kernType 'KernCompute']);
Kffs = fhandle(model.kern, model.outX, Xtest);

meanval = Kffs'*model.alpha;
if nargout>1,
    Kfsfs = fhandle(model.kern, Xtest, Xtest);
    varval = Kfsfs - Kffs'*model.Kyyinv*Kffs;
end

for d=1:model.nout, %for each output
    ymean{d} = 0;
    yvar{d} = 0;
    Ainv = mat2cell(model.Ainv, model.sizeXu, model.sizeXu);
    for q = 1:model.nlf,
        for k = [1:q-1,q+1:model.nlf],
            T = Ainv{q,k}*Kfu{d,k}.';
        end
        ymean{d} = ymean{d} + Kfu{d,q}*model.Kuuinv{q}*uast{q};
        yvar{d} = yvar{d} + ( diag(Kff{d,q})...
            - Kfu{d,q}*((model.Kuuinv{q}-Ainv{q,q})*Kfu{d,q}.' - T) );
    end
    yvar{d} = diag(yvar{d})+1/model.beta(d);
    if isfield(model,'scale') && isfield(model,'bias'),
        ymean{d} = ymean{d}*model.scale(d) + model.bias(d);
        yvar{d} = yvar{d}*model.scale(d)^2;
    end
end