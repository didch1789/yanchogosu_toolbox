function [fulldir, folder, file] = sort_ycgosu(filedir, varargin)
% sort files according to length of folder or file in input 'filedir'
% varargin(sorting criteria):
%   'folder' = sort according to folder length 
%   'file' = sort according to file length
%   'char' = out put gives a char array
% default does file sort
% output gives full directory of files, folder of files, or only files
% inspired by os.path.walk in python - !

listing = dir(filedir);

if isempty(listing)
    error('No matches for your requirement, check your input')
end

dofilesort = 1;
dofoldersort = 0;
dochar = 0;

if ~isempty(varargin)
    for i = 1:numel(varargin)
        switch varargin{i}
            case {'file'}
               dofilesort = 1;
               dofoldersort = 0; 
            case {'folder'}
               dofilesort = 0;
               dofoldersort = 1;
            case {'char'}
               dochar = 1;
            otherwise 
                error('Unidentified sorting criteria')
        end
    end
end

filelist_temp = {listing.name}'; dirlist_temp = {listing.folder}';
idx_template = 1:length(listing);

if dofilesort
    length_idx_temp = cellfun(@length, filelist_temp);
    [~, length_idx]= sort(length_idx_temp);
    idx = idx_template(length_idx);
elseif dofoldersort
    length_idx_temp = cellfun(@length, dirlist_temp);
    [~, length_idx]= sort(length_idx_temp);
    idx = idx_template(length_idx);
end

% preallocation
fulldir = cell(length(listing), 1);
folder = cell(length(listing), 1);
file = cell(length(listing), 1);

for i = 1:length(listing)
    fulldir{i, :} = fullfile(dirlist_temp{idx(i)}, filelist_temp{idx(i)});
    folder{i, :} = dirlist_temp{idx(i)};
    file{i, :} = filelist_temp{idx(i)};
end 

filt_idx = ~contains(fulldir, {'/.'});
fulldir = fulldir(filt_idx, :);
folder = folder(filt_idx, :);
file = file(filt_idx, :);

if dochar
    fulldir = char(fulldir);
    folder = char(folder);
    file = char(file);
end

end

