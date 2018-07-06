function [dKff, dKfu, dKuu, dKSigma] = dtcmgpLowerBoundGradCovMat(model)

% IBPMULTIGPLOWERBOUNDGRADCOVMAT

% IBPMULTIGP

%COPYRIGHT : Cristian Guarnizo, 2014

dKfu = cell(model.nout);
dKff = cell(model.nout,model.nlf);
dKSigma = zeros(1,model.nout);

Ainvm = model.Ainv*model.m2; %Q.Mx1
C = Ainvm*Ainvm.' + model.Ainv;
Kuuinv = blkdiag(model.Kuuinv{:});
CmKuuinv = C - Kuuinv;
for d = 1:model.nout,
    if isfield(model, 'UseMeanConstants') && model.UseMeanConstants,
        m = model.m{d} - model.mu(d);
    else
        m = model.m{d};
    end
    
    Kfdu = cell2mat( model.Kfu(d,:) );
    
    dKfu{d} = model.beta(d)*((Ainvm*m.') - (CmKuuinv)*Kfdu.').';
    
    dKff(d,:) = cellfun( @(x) -.5*x*ones(model.sizeX(d),1), num2cell(ones(1,model.nlf)*model.beta(d)), 'UniformOutput', false );
    
    Kfduqm = cell2mat(model.Psi1(d,:)');
    dKSigma(d) = sum(sum(Ainvm.* sum(Kfduqm,2) )) + .5*( -sum(sum((CmKuuinv).*model.Psi2{d}))...
        - sum(m.^2)...
        - sum(sum(cell2mat( model.Kff(d,:))))...
        + model.sizeX(d)/model.beta(d) );
end
dKfu = mat2cell(cell2mat(dKfu), cellfun('length',model.m), cellfun(@(x) size(x,1), model.latX));
dKuu = .5*(- CmKuuinv - (Kuuinv*model.P)*Kuuinv);
dKuu = mat2cell(dKuu, cellfun(@(x) size(x,1), model.latX), cellfun(@(x) size(x,1), model.latX));
induu = 1:model.nlf+1:model.nlf*model.nlf;
dKuu = dKuu(induu); %Extract only the block diagonal terms