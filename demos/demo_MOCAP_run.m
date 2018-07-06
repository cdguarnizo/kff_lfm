% DEMTOYIBPLFM Variational LFM with IBP prior over latent forces

% MULTIGP

clc
clear
close all
format short e

addpath('../dtcmgp',genpath('../globalkern'),genpath('../toolbox'),genpath('../kffkern'))

% load Golf swing data
fd = amc_to_matrix('../datasets/CMUmocap/02_03.amc');

% Remove the following channels
IndDel = [3,18,21,24,25,26,31,32,33,34,35,37,38,46];
fd(:,IndDel) = [];

outs = [2,18];
nout = length(outs);
[N, D] = size(fd);
t = (1:N)'/40; %Time stamp in seconds
scalet = 3;
%Downsample
test_ind{1} = t >= 1.7 & t <= 2.7;
test_ind{2} = t >= .7 & t <= 1.7;

y = cell(D,1);
x = cell(D,1);
xT = cell(D,1);
yT = cell(D,1);
for d = 1:D
    y{d} = fd(:, d);
    x{d} = t;
    yT{d} = y{d};
    xT{d} = x{d};
    if any(d == outs)
        ind = find(outs==d);
        y{d}(test_ind{ind}) = [];
        x{d}(test_ind{ind}) = [];
    end
end

clear fd t

% Set model Options 
options = dtcmgpOptions('dtcvar');
options.kernType = 'kfflfm';
options.optimiser = 'scg2';

options.nlf = 6;
options.numActive = 25;
options.alpha = 1;
options.NI = 200;
options.NIO = 20;
options.DispOpt = 1;
options.beta = 1e-2;
options.gamma = exp(-2);
options.isVarS = false;
options.kern.isVarS = false;
options.kern.isArd = false;
options.kern.S = 50;
options.nout = size(y,1);
ndim = size(x{1},2);
options.initialInducingPositionMethod = 'espacedInRange';

for d = 1:D
    options.bias(d) = yT{d}(1);
    options.scale(d) = std(yT{d});
end

if strcmp(options.kernType(1:2),'kf')
    name = strcat('temp/MOCAP_run_',options.kernType,'S',num2str(options.kern.S),'.mat'); 
else
    name = strcat('temp/MOCAP_run_',options.kernType,'.mat'); 
end

%% Training LFM using DTC variational approach
rng(1e5)
model = dtcmgpCreate(x, y, options);
model = dtcmgpOptimise(model, 2, 500);
params = dtcmgpExtractParam(model);
if strcmp(options.kernType(1:2),'kf')
    S = model.kern.S;
    save(name,'params','S');
else
    save(name,'params');
end
%% Plot results from the optimised model
load(name);

model = dtcmgpCreate(x, y, options);
model = dtcmgpExpandParam(model, params);
if strcmp(options.kernType(1:2),'kf')
    model.kern.S = S;
    model.kern.template.latent.S = S;
    model.kern.template.output.S = S;
end
model = dtcmgpComputeKernels(model);

[ymean, yvar] = dtcmgpPosterior(model, xT);
nmse = zeros(1,nout);
nlpd = zeros(1,nout);
fprintf('Best solution performance per output\n');
for k = 1:nout
    d = outs(k);
    ytest = yT{d}(test_ind{k});
    xtest = xT{d}(test_ind{k});
    
    nmse(k) = mysmse(ytest,ymean{d}(test_ind{k}));
    nlpd(k) = mynlpd(ytest,ymean{d}(test_ind{k}),yvar{d}(test_ind{k}));
    
    fprintf('Output %d, NMSE: %f, NLPD: %f\n',outs(k),nmse(k),nlpd(k));
    
    %Plot outputs
    figure(k)
    plotGP(ymean{d}, yvar{d}, xT{d}/scalet, y{d}, x{d}/scalet, ytest, xtest/scalet);
    grid on;
end