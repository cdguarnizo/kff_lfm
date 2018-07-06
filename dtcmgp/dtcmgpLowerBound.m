function f = dtcmgpLowerBound(model)

% IBPMULTIGPLOWERBOUND

% IBPMULTIGP
f = 0;

%% Lowerbound terms realted to u and data

% Lower bound terms related to u
Lainvm = model.Lainv'*model.m2;
f = f + .5*( (Lainvm'*Lainvm) - model.logDetA + sum(model.logDetKuu)...
    - sum( sum( model.cdq )) );

%% Lower Bound terms related to data
% Add
f = f - 0.5*sum(model.sizeX)*log(2*pi);
% Add 0.5 log Sigma_w
f = f + 0.5*sum(model.sizeX.*log(model.beta));
% Add trace (Sigma_w yy^{\top})
for d = 1:model.nout,
    f = f - 0.5*model.beta(d)*sum(model.m{d}.^2);
end

if ~isreal(f),
   warning('Imaginary part in lowerbound');
end

if isnan(f),
   warning('NaN in LB computation');
end