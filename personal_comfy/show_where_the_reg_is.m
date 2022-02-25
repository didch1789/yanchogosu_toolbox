function outputobj = find_where_the_reg_is(inputmap, varargin)

if ismac
    basedir = '/Users/jungwoo/Dropbox/Mfmri_sync';
else
    basedir = '/home/jungwoo/Dropbox/Mfmri_sync';
end

inputobj = fmri_data_rhesus(which('D99_atlas_1.2b.nii'), inputmap);
outputobj = inputobj;
outs = load(fullfile(basedir, 'D99', 'labletablesorted.mat'));

for i = 1:numel(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'regnum'
                reginfo = varargin{i+1};
            case 'regname'
                reginfo = varargin{i+1};
        end
    end
end

if ~exist('reginfo', 'var')
    error('Region number of Name should be specified!');
end

if isnumeric(reginfo)
    outputobj.dat = inputobj.dat .* (inputobj.dat == reginfo);
elseif ischar(reginfo)
    labelcell = squeeze(struct2cell(outs.labelmat));
    regNum = cell2mat(labelcell(1, contains(labelcell(2, :), reginfo)));
    outputobj.dat = inputobj.dat .* (inputobj.dat == regNum);
end




end



