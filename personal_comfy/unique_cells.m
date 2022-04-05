function [uni_cellidx] = unique_cells(cells)

k = 1;
while k <= numel(cells)
    val = cells{k};
    for i = 1:numel(cells)
        idx(i, k) = isequal(val, cells{i});
    end
    k = k + 1;
end

[~, ~, uni_cellidx] = unique(idx, 'rows', 'stable');


end