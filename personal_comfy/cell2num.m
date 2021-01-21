function out = cell2num(incell)

for i = 1:numel(incell)
    out(i, :) = str2double(incell{i});
end

end