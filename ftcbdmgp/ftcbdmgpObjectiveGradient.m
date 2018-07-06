function [f,g] = ftcbdmgpObjectiveGradient(model)


for d = 1:model.nout,
    model.Psi2{d} = 0.;
    inddq = (d-1)*model.k+1:d*model.k;
    md = model.m{d};
    Xd = model.outX{d};
    kern = model.kern;
    latX = model.latX;
    
    indXu = model.indXu;
    sizeXd = model.sizeX(d);
    BS = model.BatchSize;
    Psi2 = 0.;
    Psi1 = zeros(model.k,model.nlf);
    sKffq = zeros(1,model.nlf);
    parfor k = 1:ceil(sizeXd/BS),
        n = (k-1)*BS + 1;
        nend = n + BS - 1;
        if nend > sizeXd,
            nend = sizeXd;
        end
        indn = n:nend;
        X = Xd(indn,:);
        m = md(indn);
        
        Psi1t = zeros(size(latX{1},1),kern.nlf);
        sKffqt = zeros(1,kern.nlf);
        Kfu = zeros(length(m),length(indXu)); %n:nendxQ.M
        for q = 1:kern.nlf,
            %Compute Kff(q) and Kfu
            Kfu = globalKfuCompute(kern, X, latX{q}, d, q);
            Kffq = globalKffCompute(kern, X, d, q);
            
            Psi1t(:,q) = Psi1t(:,q) + Kfu'*m; %Mx1
            sKffqt(q) = sKffqt(q) + sum(Kffq);
            Kfut(:,indXu == q) = Kfu;
        end
        Psi2t = Kfut'*Kfut; %Q.MxQ.M
        
        sKffq = sKffq + sKffqt;
        Psi1 = Psi1 + Psi1t;
        Psi2 = Psi2 + Psi2t; %Q.MxQ.M
    end
    model.sKffq(d,:) = sKffq;
    model.Psi1(inddq,:) = Psi1;
    model.Psi2{d} = Psi2;
end

Lainvm2 = Lainv'*model.m2;
% LowerBound computation
f = -.5*(sum(model.beta.*model.m2d) + Lainvm2'*Lainvm2 - model.logDetA + sum(model.logDetKuu)... 
    - sum(sum(model.cdq.*model.etadq)) );

%% Phase 2 - Compute the gradients
Ainvm = model.Ainv*model.m2; %Q.Mx1
C = Ainvm*Ainvm.' + model.Ainv;
Kuuinv = blkdiag(model.Kuuinv{:});
CmKuuinv = C - Kuuinv;

gKern = zeros(1,model.kern.nParams);
gBeta = zeros(1,model.nout);
for d = 1:model.nout,
    indd = d*ones(1,model.sizeX(d));
    EZ2 = model.etadq(d,:)'*model.etadq(d,:) - diag(model.etadq(d,:).^2) + diag(model.etadq(d,:));
    EZ22 = EZ2(model.indXu,model.indXu);
    EZ = model.etadq(indd,model.indXu);
    
    C1 = model.beta(d)*(EZ'.*(Ainvm*m'));
    C2p = EZ22.*CmKuuinv;
    C2 = model.beta(d)*C2p;
    
    n = 1;
    while n < model.sizeX(d),
        nend = n + model.BatchSize;
        if nend > model.sizeX(d),
            nend = model.sizeX(d);
        end
        indn = n:nend;
        X = model.outX{d}(indn,:);
        m = model.m{d}(indn);
        
        for q = 1:model.nlf,
            %Compute Kfu
            fhandle = str2func([kern.type 'KfuCompute']);
            Kfu = fhandle(model.kern,X,model.latX{q},d,q);
            
            dKfu = (C1 - C2*Kfu.').';
            
            dKff = -.5*model.etadq(d,q)*model.beta(d)*eye(length(m));
            
            fhandle = str2func([kern.type 'KffKfuGradient']);
            gKern = gKern + fhandle(model.kern,X,model.latX{q},dKff,dKfu,d,q);
        end
    end
    
    %Beta gradients
    Kfduqm = cell2mat(model.Psi1(d,:)');
    gBeta(d) = sum(sum(Ainvm.* sum(model.etadq(d,model.indXu)'.*Kfduqm,2) ))...
        +.5*( -sum(sum(C2p.*model.Psi2{d}))...
        - model.md2(d)+ model.sizeX(d)/model.beta(d) - sum(model.etadq(d,:)*sKffq(d,:)));
end

factors = kernFactors(model.kern, 'gradfact');
gKern(factors.index) = gKern(factors.index).*factors.val;

for q = 1:model.nlf,
    dKuu = -.5*(CmKuuinv(indq,indq) + model.Kuuinv{q}*model.P(indq,indq)*model.Kuuinv{q});
    fhandle = str2func([kern.type 'KuuGradient']);
    gKern = gKern + fhandle(model.kern, model.latX{q}, dKuu, q);
end

fhandle = str2func([model.betaTransform 'Transform']);
gBeta = gBeta.*fhandle(model.beta, 'gradfact');

g = [gKern,gGamma,gBeta];
