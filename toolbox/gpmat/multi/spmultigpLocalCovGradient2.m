function [dLdKyy, dLdKuy, dLdKuu, dLdmu, dLdbeta] = spmultigpLocalCovGradient2(model)

% SPMULTIGPLOCALCOVGRADIENT Computes the derivatives of the likelihood
% objective function corresponding to the sparse approximation, with
% respect to the kernels of the multi output gp.
%
% FORMAT
% DESC Computes the derivatives of the likelihood objective function corresponding
% to the sparse approximation
% ARG model : the model structure containing the information about
% the model.
% RETURN dLdKyy: the derivatives of the likelihood with respect to the
% kernels Kyy.
% RETURN dLdKuy: the derivatives of the likelihood with respect to the
% kernels Kuy.
% RETURN dLdKuu: the derivatives of the likelihood with respect to the
% kernels Kuu.
% RETURN dLddmu: the derivatives of the likelihood with respect to the
% mean function.
% RETURN dLdbeta: the derivatives of the likelihood with respect to the
% noise beta parameters.

% COPYRIGHT : Mauricio A Alvarez, 2008, 2009

% MULTIGP

Ainv = cell2mat(model.Ainv);
AinvKuyDinvy = cell2mat(model.AinvKuyDinvy);
Kyu = cell2mat(model.Kyu);

% /~MAURICIO : Last code
% C = mat2cell(Ainv + (AinvKuyDinvy*AinvKuyDinvy'), model.k*ones(1,model.nlf), model.k*ones(1,model.nlf));
% CKuy = mat2cell(cell2mat(C)*Kyu', model.k*ones(1,model.nlf), cellfun('length',model.m));
% ~/

CC = Ainv + (AinvKuyDinvy*AinvKuyDinvy');
C = mat2cell(CC, model.k, model.k);
[sqrtCC, jitter] = jitChol(CC);
sqrtC = mat2cell(sqrtCC, model.k, model.k);
CKuy = mat2cell(CC*Kyu', model.k, cellfun('length',model.m));
sqrtCKuy = mat2cell(sqrtCC*Kyu', model.k, cellfun('length',model.m));

switch model.approx
    case 'dtc'
        dLdKyy = cell(model.nout,1);
        dLdKuy = cell(model.nlf,model.nout);
        for r =1:model.nlf,
            for k= 1:model.nout,
                if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt)
                    switch model.noiseOpt
                        case 0
                            dLdKuy{r,k} = - CKuy{r,k}*model.beta(k) + model.AinvKuyDinvy{r}*model.m{k}'*model.beta(k);
                        case 1                            
                        if isfield(model, 'nRepeats') && ~isempty(model.nRepeats)
                            dLdKuy{r,k} =  - CKuy{r,k}*diag(model.beta(k)*model.nRepeats{k}) ...
                                + model.AinvKuyDinvy{r}*model.m{k}'*diag(model.beta(k)*model.nRepeats{k});                            
                        end
                        case 2
                        case 3
                            dLdKuy{r,k} =  - CKuy{r,k}*diag(model.beta{k})...
                                + model.AinvKuyDinvy{r}*model.m{k}'*diag(model.beta{k});                            
                    end
                else
                    dLdKuy{r,k} = - CKuy{r,k}*model.beta(k) + model.AinvKuyDinvy{r}*model.m{k}'*model.beta(k);
                end
            end
        end
        dLdKuu = cell(model.nlf,1);
        for r =1:model.nlf,
            dLdKuu{r} = 0.5*((model.Kuuinv{r} - C{r,r}));
        end
        dLdbeta = zeros(1,model.nout);
        dLdmu = cell(model.nout,1);
        for k =1:model.nout,
            DinvKyuAinvKuyDinvy = zeros(size(model.X{k+model.nlf},1),1);
            H = zeros(size(model.X{k+model.nlf},1));
            for r =1:model.nlf,
                H = H + model.Kyu{k,r}*model.AinvKuyDinvy{r}*model.m{k}'+ (model.Kyu{k,r}*model.AinvKuyDinvy{r}*model.m{k}')'...
                    - model.Kyu{k,r}*CKuy{r,k};
                DinvKyuAinvKuyDinvy = DinvKyuAinvKuyDinvy + model.KuyDinv{r,k}'*model.AinvKuyDinvy{r};
            end            
            if ~strcmp(model.kernType,'gg')
                if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt)
                    switch model.noiseOpt
                        case 0
                            dLdmu{k} = model.beta(k)*model.m{k} - DinvKyuAinvKuyDinvy;
                        case 1
                            dLdmu{k} = (model.beta(k)*model.nRepeats{k}).*model.m{k} - DinvKyuAinvKuyDinvy;
                        case 2
                        case 3
                            dLdmu{k} = (model.beta{k}).*model.m{k} - DinvKyuAinvKuyDinvy;
                    end                    
                else
                    dLdmu{k} = model.beta(k)*model.m{k} - DinvKyuAinvKuyDinvy;
                end
            end
            if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt)
                switch model.noiseOpt
                    case 0
                        dLdbeta(k) = -0.5*( - (size(model.X{k+model.nlf},1)*1/model.beta(k) - ...
                            model.m{k}'*model.m{k} + trace(H)));
                    case 1
                        if isfield(model, 'nRepeats') && ~isempty(model.nRepeats)
                            dLdbeta(k) = -0.5*model.beta(k)^(-2)*trace(sparseDiag((1./model.nRepeats{k}))*...
                                (model.Dinv{k}*( - ...
                                (model.D{k} - model.m{k}*model.m{k}' + H))*model.Dinv{k}));
                        end
                    case 2
                    case 3
                        dLdbeta = [];
                end
            else
                dLdbeta(k) = -0.5*(- (size(model.X{k+model.nlf},1)*1/model.beta(k)...
                    - model.m{k}'*model.m{k} + trace(H)));
            end
        end
    case {'fitc','pitc'}
        Q = cell(model.nout,1);
        sqrtQ = cell(model.nout,1);
        dLdKyy = cell(model.nout,1);
        dLdmu = cell(model.nout,1);
        for k =1:model.nout,
            DinvKyuAinvKuyDinvy = zeros(size(model.X{k+model.nlf},1),1);
            Q{k} = zeros(size(model.X{k+model.nlf},1));
            for r =1:model.nlf,
                Q{k} = Q{k} + model.Kyu{k,r}*model.AinvKuyDinvy{r}*model.m{k}'...
                    + (model.Kyu{k,r}*model.AinvKuyDinvy{r}*model.m{k}')' - (sqrtCKuy{r,k})'*sqrtCKuy{r,k};
                DinvKyuAinvKuyDinvy = DinvKyuAinvKuyDinvy + model.KuyDinv{r,k}'*model.AinvKuyDinvy{r};
            end
            Q{k} = model.Dinv{k}*(model.D{k} - model.m{k}*model.m{k}' + Q{k})*model.Dinv{k};
            if strcmp(model.approx, 'fitc'),
                Q{k} = sparseDiag(diag(Q{k}));
            end
            %[sqrtQ{k}, jitter] = jitChol(Q{k});
            dLdKyy{k} = -0.5*Q{k};
            if ~strcmp(model.kernType,'gg')
                dLdmu{k} = model.Dinv{k}*model.m{k} - DinvKyuAinvKuyDinvy;
            end
        end
        dLdKuy = cell(model.nlf,model.nout);
        for r =1:model.nlf,
            for k= 1:model.nout,
                dLdKuy{r,k} = model.KuuinvKuy{r,k}*Q{k}...
                    - CKuy{r,k}*model.Dinv{k} + model.AinvKuyDinvy{r}*model.m{k}'*model.Dinv{k};
            end
        end
        dLdKuu = cell(model.nlf,1);
        for r =1:model.nlf,
            KuuinvKuyQKyuKuuinv = zeros(model.k(r));
            dLdKuu{r} = zeros(model.k(r));
            for k=1: model.nout,
                %KuuinvKuyQKyuKuuinv =  KuuinvKuyQKyuKuuinv + model.KuuinvKuy{r,k}*Q{k}*model.KuuinvKuy{r,k}';
                KuuinvKuyQKyuKuuinv =  KuuinvKuyQKyuKuuinv + model.KuuinvKuy{r,k}*Q{k}*model.KuuinvKuy{r,k}';
            end
            dLdKuu{r} = 0.5*((model.Kuuinv{r} - C{r,r}) - KuuinvKuyQKyuKuuinv);
        end
        if nargout>4
            dLdbeta = zeros(1,model.nout);
            for k =1:model.nout,
                if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt)
                    switch model.noiseOpt
                        case 0
                            dLdbeta(k) = 0.5*model.beta(k)^(-2)*trace(Q{k});
                        case 1
                            if isfield(model, 'nRepeats') && ~isempty(model.nRepeats)
                                dLdbeta(k) = 0.5*model.beta(k)^(-2)*trace(diag((1./model.nRepeats{k}))*Q{k});
                            end
                        case 2
                            %
                        case 3
                            dLdbeta = [];
                    end
                else
                    dLdbeta(k) = 0.5*model.beta(k)^(-2)*trace(Q{k});
                end
            end
        end
    case 'dtcvar'
        dLdKyy = cell(model.nout,1);
        for k=1:model.nout,
            dLdKyy{k} =  -0.5*model.Dinv{k};
        end
        dLdKuy = cell(model.nlf,model.nout);
        for r =1:model.nlf,
            for k= 1:model.nout,
                if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt)
                    switch model.noiseOpt
                        case 0
                            dLdKuy{r,k} = model.KuuinvKuy{r,k}*model.beta(k)...
                                - CKuy{r,k}*model.beta(k) + model.AinvKuyDinvy{r}*model.m{k}'*model.beta(k);
                        case 1
                            if isfield(model, 'nRepeats') && ~isempty(model.nRepeats)
                                dLdKuy{r,k} = model.KuuinvKuy{r,k}*diag(model.beta(k)*model.nRepeats{k})...
                                    - CKuy{r,k}*diag(model.beta(k)*model.nRepeats{k}) ...
                                    + model.AinvKuyDinvy{r}*model.m{k}'*diag(model.beta(k)*model.nRepeats{k});

                            end
                        case 2
                            %
                        case 3
                            dLdKuy{r,k} = model.KuuinvKuy{r,k}*diag(model.beta{k}) - CKuy{r,k}*diag(model.beta{k})...
                                + model.AinvKuyDinvy{r}*model.m{k}'*diag(model.beta{k});
                    end
                else
                    dLdKuy{r,k} = model.KuuinvKuy{r,k}*model.beta(k)...
                        - CKuy{r,k}*model.beta(k) + model.AinvKuyDinvy{r}*model.m{k}'*model.beta(k);
                end
            end
        end
        dLdKuu = cell(model.nlf,1);
        for r =1:model.nlf,
            KuuinvKuyQKyuKuuinv = zeros(model.k(r));
            dLdKuu{r} = zeros(model.k(r));
            for k=1: model.nout,
                KuuinvKuyQKyuKuuinv =  KuuinvKuyQKyuKuuinv + model.KuuinvKuy{r,k}*model.Dinv{k}*model.KuuinvKuy{r,k}';
            end
            dLdKuu{r} = 0.5*((model.Kuuinv{r} - C{r,r}) - KuuinvKuyQKyuKuuinv);
        end
        dLdbeta = zeros(1,model.nout);
        dLdmu = cell(model.nout,1);
        for k =1:model.nout,
            DinvKyuAinvKuyDinvy = zeros(size(model.X{k+model.nlf},1),1);
            H = zeros(size(model.X{k+model.nlf},1));
            for r =1:model.nlf,
                H = H + model.Kyu{k,r}*model.AinvKuyDinvy{r}*model.m{k}'+ (model.Kyu{k,r}*model.AinvKuyDinvy{r}*model.m{k}')'...
                    - model.Kyu{k,r}*CKuy{r,k};
                DinvKyuAinvKuyDinvy = DinvKyuAinvKuyDinvy + model.KuyDinv{r,k}'*model.AinvKuyDinvy{r};
            end
            if ~strcmp(model.kernType,'gg')
                if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt)
                    switch model.noiseOpt
                        case 0
                            dLdmu{k} = model.beta(k)*model.m{k} - DinvKyuAinvKuyDinvy;
                        case 1
                            dLdmu{k} = (model.beta(k)*model.nRepeats{k}).*model.m{k} - DinvKyuAinvKuyDinvy;
                        case 2
                        case 3
                            dLdmu{k} = (model.beta{k}).*model.m{k} - DinvKyuAinvKuyDinvy;
                    end
                else
                    dLdmu{k} = model.beta(k)*model.m{k} - DinvKyuAinvKuyDinvy;
                end
            end
            if isfield(model, 'noiseOpt') && ~isempty(model.noiseOpt)
                switch model.noiseOpt
                    case 0
                        dLdbeta(k) = -0.5*(sum(model.Ktilde{k}) - (size(model.X{k+model.nlf},1)*1/model.beta(k) - ...
                            model.m{k}'*model.m{k} + trace(H)));
                    case 1
                        if isfield(model, 'nRepeats') && ~isempty(model.nRepeats)
                            dLdbeta(k) = -0.5*model.beta(k)^(-2)*trace(sparseDiag((1./model.nRepeats{k}))*...
                                (model.Dinv{k}*(sparseDiag(model.Ktilde{k}) - ...
                                (model.D{k} - model.m{k}*model.m{k}' + H))*model.Dinv{k}));
                        end
                    case 2
                    case 3
                        dLdbeta = [];
                end
            else
                dLdbeta(k) = -0.5*(sum(model.Ktilde{k}) - (size(model.X{k+model.nlf},1)*...
                    1/model.beta(k) - model.m{k}'*model.m{k} + trace(H))  );
            end
        end
end