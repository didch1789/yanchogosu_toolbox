function mkdir_BIDS_format
% Run this function in "~/data/$prjname"
basedir = pwd;

mkdir(fullfile(basedir, 'audiovideo'))
mkdir(fullfile(basedir, 'behavioral'))
mkdir(fullfile(basedir, 'documents'))
mkdir(fullfile(basedir, 'imaging'))
mkdir(fullfile(basedir, 'imaging/preprocessed'))
mkdir(fullfile(basedir, 'imaging/raw'))
mkdir(fullfile(basedir, 'imaging/dicom_from_scanner'))
mkdir(fullfile(basedir, 'physio'))
mkdir(fullfile(basedir, 'figures'))

end