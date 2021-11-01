function r_squared = calc_rsquared_ycgosu(X, y, varargin)
% calculate r_squared based on ordinary least squared method.
%
% Inputs
%   X   : matrix with #observations X #features
%   y   : matrix with #observations X 1
%   varargin:
%       'adjusted' : calculate adjusted r_squared.

doadjusted = false;

for i = 1:numel(varargin)
    if isa(varargin{i}, 'char') 
        switch varargin{i}
            case 'adjusted'
                doadjusted = true;
        end
    end
end

[N, p] = size(X);

b = glmfit(X, y);
yfit = b(1) + X * b(2:end);


if doadjusted
    orir_squared = 1 - (((y - yfit)' * (y - yfit)) / ((y-mean(y))' * y-mean(y)));
    r_squared = 1 - ((1 - orir_squared) * (N - 1) / (N - p - 1));
else
    r_squared = 1 - (((y - yfit)' * (y - yfit)) / ((y-mean(y))' * y-mean(y)));
end

end