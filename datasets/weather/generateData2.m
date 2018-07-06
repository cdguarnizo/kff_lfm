% This code builds the data for weather data
load bra
ybra = y(:,4);
load cam
ycam = y(:,4);
load chi
ychi = y(:,4);
load sot
ysot = y(:,4);

range = x >= 0 & x <= 15;
xtemp = x(range);

ytemp = [ybra, ycam, ychi, ysot];
ytemp = ytemp(range, :);
clear ybra ycam ychi ysto;
% ymean = mean(y);
% ystd = std(y);
% ytemp = ytemp - repmat(ymean, size(ytemp,1), 1);
% ytemp = ytemp ./ repmat(ystd, size(ytemp,1), 1);

ytest = cell(4,1);
y = ytest;
Xtest = cell(4,1);
X = Xtest;
for d=1:size(ytemp,2),
    inddel = ytemp(:,d) ~= -1;
    ytest{d} = ytemp(inddel,d);
    Xtest{d} = xtemp(inddel);
    if d == 2,
        range = Xtest{d} < 6.25 | Xtest{d} > 8.25;
        y{d} = ytest{d}(range);
        X{d} = Xtest{d}(range);
    elseif d == 3,
        range = Xtest{d} < 13.5 | Xtest{d} > 14.2;
        y{d} = ytest{d}(range);
        X{d} =  Xtest{d}(range);
    else
        X{d} = Xtest{d};
        y{d} = ytest{d};
    end
end

x = X;
xT = Xtest;
yT = ytest;
save weatherdatafull x y xT yT