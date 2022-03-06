function binned = bin_raster(rasterdat, binsize)
% Inputs
%   rasterdat : raw spike data in 1ms. (neurons X timeseries)
%   binsize   : time bin size in ms (e.g. 10) 
%
% Outputs
%   binned: binned (summed wihthin bin) raster 
% e.g. 
%   size(rasterdat) = 100 200 (100 neurons' 200 ms data)
%   binsize = 10;
%   size(binned)    = 100 20
%
%   size(rasterdat) = 100 205 (100 neurons' 200 ms data)
%   binsize = 10
%   size(binned)    = 100 20 (first 5 bin will contain 11 binsize)
neuronnum   = size(rasterdat, 1); 
timesize    = size(rasterdat, 2);
Q           = fix(timesize / binsize);
R           = mod(timesize, binsize);
binned      = NaN(neuronnum, Q);

l = 1;
for bb = 1:Q
    if R ~= 0
        k = binsize + 1;
        R = R - 1;
    else
        k = binsize;
    end
    
    binned(:, bb) = sum(rasterdat(:, l:l+k-1), 2);
    l = l + k;

end