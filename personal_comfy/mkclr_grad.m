function out = mkclr_grad(index, targ_clrs)

if numel(index) ~= numel(targ_clrs)
    error('Input size of 1st and 2nd inputs are not identical');
end

clridx = mat2cell(1:sum(index), 1, index);

for ii = 1:numel(clridx)
    out(clridx{ii}, :) = clr_interpolate(targ_clrs{ii}(1, :), targ_clrs{ii}(2, :), numel(clridx{ii}));
end

end