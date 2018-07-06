function f = dtcmgpLowerBoundOptimise(model)

% IBPMULTIGPLOWERBOUND FOR IBP Parameters

% IBPMULTIGP
f = 0;

%% Lowerbound terms realted to u and data

% Lower bound terms related to u
f = f + .5*( (model.m2'*(model.A\model.m2)) - model.logDetA + sum(model.logDetKuu)...
    - sum(sum( model.cdq.*model.etadq )) );
        
%% Lower Bound terms related to data
% Add 0.5 log Sigma_w
f = f + 0.5*sum(model.sizeX.*log(model.beta));
% Add trace (Sigma_w yy^{\top})
f = f - 0.5*sum(model.beta(model.indX)'.*cell2mat(model.m).^2);