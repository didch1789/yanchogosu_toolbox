function idx = depcols_ycgosu(X, varargin)
% varagin:
%    tol: tolerance default tol is 1E-5
tol = 1E-5;

for i = 1:numel(varargin)
    switch varargin{i}
        case 'tol'
            tol = varargin{i+1};
    end
end

issame = zeros(size(X, 2), size(X, 2)-1); 

for i = 1:size(X, 2)
    issame(i, :) = corr(X(:, i),X(:, setdiff(1:size(X, 2), i)));
end

[a, b] = find(issame == tol);

idx = [a(1:end/2, :) b(1:end/2, :)];

end