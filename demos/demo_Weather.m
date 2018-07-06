% Demo using DTCVAR for learning Weather Data

% MULTIGP

clc
clear
close all
format short e

addpath('../dtcmgp',genpath('../globalkern'),genpath('../toolbox'),genpath('../kffkern'))
% load data
load ../datasets/weather/weatherdata.mat
D = size(y,1);

names = {'Bramblemet','Cambermet','Chimet','Sotonmet'};
test_ind{1} = xT{2} >= 10.2 & xT{2} <= 10.8;
test_ind{2} = xT{3} >= 13.5 & xT{3} <= 14.2;
outs = [2,3];
nout = length(outs);

% Set the Options 
options = dtcmgpOptions('dtcvar');
options.kernType = 'kffsim';
options.optimiser = 'scg2';

options.nlf = 6;
options.numActive = 200;
options.alpha = 1;
options.NI = 200;
options.NIO = 20;
options.DispOpt = 1;
options.beta = 1e-2;
options.gamma = exp(-2);
options.isVarS = false;
options.kern.isVarS = false;
options.kern.isArd = false;
options.kern.S = 10;
options.nout = size(y,1);
ndim = size(x{1},2);
options.initialInducingPositionMethod = 'espacedInRange';

if strcmp(options.kernType(1:2),'kf'),
    name = strcat('temp/Weather',options.kernType,'S',num2str(options.kern.S),'.mat'); 
else
    name = strcat('temp/Weather',options.kernType,'.mat'); 
end

%% Training IBPLFM variational approach
rng(1e6)
model = dtcmgpCreate(x, y, options);
model = dtcmgpOptimise(model, 2, 500);
params = dtcmgpExtractParam(model);
if strcmp(options.kernType(1:2),'kf'),
    S = model.kern.S;
    save(name,'params','S');
else
    save(name,'params');
end
%% Plot results from the optimised model
load(name);

model = dtcmgpCreate(x, y, options);
model = dtcmgpExpandParam(model, params);
if strcmp(options.kernType(1:2),'kf'),
    model.kern.S = S;
    model.kern.template.latent.S = S;
    model.kern.template.output.S = S;
end
model = dtcmgpComputeKernels(model);

[ymean, yvar] = dtcmgpPosterior(model, xT);

nmse = zeros(1,nout);
nlpd = zeros(1,nout);
fprintf('Best solution performance per output\n');
for k = 1:nout,
    d = outs(k);
    ytest = yT{d}(test_ind{k});
    xtest = xT{d}(test_ind{k});
    
    nmse(k) = mysmse(ytest,ymean{d}(test_ind{k}));
    nlpd(k) = mynlpd(ytest,ymean{d}(test_ind{k}),yvar{d}(test_ind{k}));
    
    fprintf('Output %d, NMSE: %f, NLPD: %f\n',outs(k),nmse(k),nlpd(k));
    
    %Plot outputs
    figure(k)
    plotGP(ymean{d}, yvar{d}, xT{d}, y{d}, x{d}, ytest, xtest);
    %title(names{d})
end