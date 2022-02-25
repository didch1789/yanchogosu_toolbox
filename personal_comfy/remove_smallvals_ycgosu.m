function obj = remove_smallvals_ycgosu(obj, threshold)
% this function delete all the values below threshold in the fmri_data
% object.
idx = ~(obj.dat < threshold);
obj2 = obj;

obj2.dat = obj2.dat(idx);
obj2.removed_voxels = false(size(obj2.dat, 1), 1);

obj2.volInfo.wh_inmask = obj2.volInfo.wh_inmask(idx, :); 
obj2.volInfo.n_inmask = sum(idx);
obj2.volInfo.xyzlist = obj2.volInfo.xyzlist(idx, :); 
obj2.volInfo.cluster = obj2.volInfo.cluster(idx, :);
numidx = find(obj2.volInfo.image_indx == 1);
numidx2 = numidx(idx);

tempidx = zeros(obj2.volInfo.nvox, 1);
tempidx(numidx2) = 1;
obj2.volInfo.image_indx = logical(tempidx);
obj2.volInfo.fname = [obj2.volInfo.fname filesep sprintf('remove vals below %d', threshold)];
obj = obj2;

end