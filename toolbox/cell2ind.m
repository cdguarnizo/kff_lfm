function Xind = cell2ind(X)
% This funcion allow to convert from cell to index scheme

Xind.val = cell2mat(X);
Xind.ind = zeros(size(X.val));
endVal = 0;
for d = 1:size(X,1),
    startVal = endVal + 1;
    endVal = starVal + size(X{d},2);
    Xind.ind(startVal:endVal) = d;
end