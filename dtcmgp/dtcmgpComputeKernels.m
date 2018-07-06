function model = dtcmgpComputeKernels(model)

% IBPMULTIGPCOMPUTEKERNELS

% IBPMULTIGP

fhandle = str2func([model.kernType 'KernCompute']);
if isfield(model, 'gamma') && ~isempty(model.gamma)
    [model.Kff, model.Kfu, model.Kuu] = fhandle(model.kern, ...
        model.outX, model.latX, model.gamma);
else
    [model.Kff, model.Kfu, model.Kuu] = fhandle(model.kern, ...
        model.outX, model.latX);
end

% Compute Kuuinv
model.m2 = cell(model.nlf,1);
for q = 1:model.nlf,
    Lu = jitChol(model.Kuu{q});
    model.Kuu{q} = Lu.'*Lu;
    model.logDetKuu(q) = 2.*sum(log(diag(Lu)));
    Lu = Lu\eye(size(Lu));
    model.Kuuinv{q} = Lu*Lu';
    
    model.m2{q} = 0.;
end

% Compute cdq and Psi
model.cdq = zeros(model.nout,model.nlf);
model.Psi2 = cell(1, model.nout);
model.Psi1 = cell(model.nout, model.nlf);
model.P = 0.;
for d = 1:model.nout,
    Kfdu = cell2mat( model.Kfu(d,:) );
    model.Psi2{d} = Kfdu'*Kfdu; %Q.MxQ.M
    for q = 1:model.nlf,
        indq = model.indXu == q;
        Psi2dq = model.Psi2{d}(indq,indq);
        model.cdq(d,q) = model.beta(d)*( sum(model.Kff{d,q}) - ...
            sum(sum( model.Kuuinv{q}.*Psi2dq )) );
        model.Psi1{d,q} = model.Kfu{d,q}'*model.m{d}; %Mx1
        model.m2{q} = model.m2{q} + model.beta(d)*model.Psi1{d,q}; %Mx1
    end
    model.P = model.P + model.beta(d)*(model.Psi2{d}); %Q.MxQ.M
end
model.m2 = cell2mat(model.m2);

model.A = blkdiag(model.Kuu{:}) + model.P;  %Q.MxQ.M
[La, jitter] = jitChol(model.A);
if jitter > model.minJit,
    error('Jitter higher than tolerance.')
end
model.A = La.'*La;
model.logDetA = 2.*sum(log(diag(La)));
model.Lainv = La\eye(size(La));
model.Ainv = model.Lainv*model.Lainv.';