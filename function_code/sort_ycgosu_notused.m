% function sorted = sort_ycgosu(filedir, varargin)

% sort files in ascending order based on filenames' length. excluding
% files that starts with '.' 

% options
% 'mat' = contains only mat files in the output
% ('idx', logical_idx) = order files according to subject index, should be
% used with 'mat'

dofile = 0;
doidx = 0;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'mat', 'nii', 'png', 'jpeg', 'dicom', 'tar'}
                dofile = 1;
                form = varargin{i};
            case {'indexing', 'idx'}
                doidx = 1;
                idx = varargin{i+1};
        end
    end
end

if ~exist(filedir, 'dir'), error('file directory does not exist.'), end
if ~isa(filedir, 'char'), error('Inadequate type. string file directory should be the input.')
    

else 
    listing = dir(filedir);
    filelist_temp = {listing.name};
    dir_idx_temp = {listing.isdir};
    dir_idx = cell2mat(dir_idx_temp);

    filelist_temp = filelist_temp(~dir_idx);

    remain_idx = ones(1, numel(filelist_temp));
    for i = 1:numel(filelist_temp)
        if startsWith(filelist_temp{i}, '.')
            remain_idx(i) = 0;
        end
    end

    filelist = filelist_temp(logical(remain_idx));
    filelist = sort(filelist);

    for i = 1:numel(filelist)
        length_vals(i, :) = length(filelist{i});
    end

    uniques = unique(length_vals);
    unique_mtrx = zeros(numel(filelist), numel(uniques));
    sorted = {};

    k = 1;
    for i = 1:numel(uniques)
        unique_mtrx(:, i) = (length_vals == uniques(i));
        temp = filelist(logical(unique_mtrx(:, i)));
        sorted(k:k+numel(temp)-1) = temp';
        k = k + numel(temp);
    end
    sorted = sorted';
    
    if dofile
        for i = 1:numel(sorted)
            if endsWith(sorted{i}, strcat('.', form))
                fileidx(i) = 1;
            end
        end
        sorted = sorted(logical(fileidx));
    end
    
    if dofile && doidx
        while sum(idx) < numel(sorted)
            disp(sorted)
            x = input('type rows to be deleted:  \n');
            sorted(x) = [];
        end
        
        sorted2 = num2cell(NaN(numel(idx), 1));
        k = 1;
        
        for i=1:numel(idx)
            if idx(i) == 1
                sorted2(i) = sorted(k);
                k = k + 1;
            end
        end
        sorted = sorted2;
    end
end
end

