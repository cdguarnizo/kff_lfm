function model = ftcmgpComputeKernel(model)

% FTCMGPCOMPUTEKERNELS
% FTCMGP

fhandle = str2func([model.kernType 'KernCompute']);
model.Kyy = fhandle(model.kern, model.outX);

%model.Kyy = zeros(size(Kffq));
% for k = 1:model.nlf,
%     SS = model.S(:, k)'*model.S(:, k);
%     model.Kyy = model.Kyy + SS(model.outX.ind, model.outX.ind).*Kffq{k};
% end

if model.includeNoise,
    % Generating noise vector
    noise = cell(model.nout,1);
    for d = 1:model.nout,
        noise{d} = zeros(model.sizeX(d),1);
        noise{d}(:) = 1/model.beta(d);
    end
    noise = cell2mat(noise);
    % Adding noise
    [row, col] = size(model.Kyy);
    diagind = 1:row+1:row*col;
    model.Kyy(diagind) = model.Kyy(diagind) + noise';
end
% Update alpha
yvec = cell2mat(model.y);
model.L = jitChol(model.Kyy);
%model.Kyyinv = model.L'\(model.L\eye(size(model.Kyy)));
Li = model.L\eye(size(model.Kyy,1));
model.Kyyinv = Li*Li'; % pdinv(model.Kyy)
%model.alpha =  model.L'\(model.L\yvec); % model.Kyyinv*yvec;
model.alpha =  model.Kyyinv*yvec;