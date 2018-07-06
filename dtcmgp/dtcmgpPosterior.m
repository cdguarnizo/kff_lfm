function [ymean, yvar] = dtcmgpPosterior(model, Xtest)

ymean = cell(model.nout, 1);
yvar = cell(model.nout, 1);

fhandle = str2func([model.kernType 'KernCompute']);
if isfield(model, 'gamma') && ~isempty(model.gamma)
    [Kff, Kfu, ~] = fhandle(model.kern, Xtest, model.latX, model.gamma);
else
    [Kff, Kfu, ~] = fhandle(model.kern, Xtest, model.latX);
end
Kuu = blkdiag(model.Kuu{:});
temp = Kuu*model.Ainv;
uast = temp*model.m2;
uast = mat2cell(uast, model.sizeXu);

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
        yvar{d} = yvar{d}*model.scale(d)*model.scale(d);
    end
end