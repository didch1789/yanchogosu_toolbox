function out = find_biggest_region(obj, varargin)

regobj = region(obj);
voxnums = zeros(numel(regobj), 1);
for ni = 1:numel(regobj)
    voxnums(ni, :) = regobj(ni).numVox;
end
[~, maxidx] = max(voxnums);
regout = regobj(maxidx);
out = region2imagevec(regout);


end