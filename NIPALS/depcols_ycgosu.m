function idx = depcols_ycgosu(X)

issame = zeros(size(X, 2), size(X, 2)-1); 

for i = 1:size(X, 2)
    issame(i, :) = corr(X(:, i),X(:, setdiff(1:size(X, 2), i)));
end

[a, b] = find(issame == 1);

idx = [a(1:end/2, :) b(1:end/2, :)];

end