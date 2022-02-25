function [r_squared, y_error] = calc_rsquared_ycgosu(X, y, varargin)
% calculate r_squared based on ordinary least squared method.
%
% Inputs
%   X   : matrix with #observations X #features
%   y   : matrix with #observations X 1
%   varargin:
%       'adjusted' : calculate adjusted r_squared.

do_adjusted = false;
do_crossval = false;
do_addIntcpt= false;

for i = 1:numel(varargin)
    if isa(varargin{i}, 'char') 
        switch varargin{i}
            case 'adjusted'
                do_adjusted = true;
            case 'nfolds'
                do_crossval = true;
                nfolds = varargin{i+1};
            case 'addIntcpt'
                do_addIntcpt = true;
        end
    end
end

[N, p] = size(X);
yfit = NaN(N, 1);
if do_crossval
    CVs = cvpartition(N, 'KFold', nfolds);
    fprintf('Doing %02d folds CV...\n', nfolds);
    for nf = 1:nfolds
        TrIdx = CVs.training(nf);   TeIdx = CVs.test(nf);
        Xtrain = X(TrIdx, :);       ytrain = y(TrIdx, :);
        Xtest  = X(TeIdx, :);       
        if do_addIntcpt
            newX = [ones(size(Xtrain, 1), 1) Xtrain];
            b = pinv(newX) * ytrain;
            yfit(TeIdx) = b(1) + Xtest * b(2:end);
        else
            b = pinv(Xtrain) * ytrain;
            yfit(TeIdx) = Xtest * b;
        end
    end
else
    if do_addIntcpt
        newX = [ones(size(X, 1), 1) X];
        b = pinv(newX) * y;
        yfit = newX * b;
    else
        b = pinv(X) * y;
        yfit = X * b;
    end
end

y_error = y - yfit;

if do_adjusted
    orir_squared = 1 - (((y - yfit)' * (y - yfit)) / ((y-mean(y))' * y-mean(y)));
    r_squared = 1 - ((1 - orir_squared) * (N - 1) / (N - p - 1));
else
    r_squared = 1 - (((y - yfit)' * (y - yfit)) / ((y-mean(y))' * y-mean(y)));
end

end