function out = perm_corr(x,y)
rng('default');
for i = 1:10000
 idx_x = randperm(numel(x));
 idx_y = randperm(numel(x));
 r(i,1) = corr(x(idx_x), y(idx_y));
end

out.permP = sum(corr(x,y) < r)./10000;

end