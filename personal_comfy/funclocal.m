function funclocal(niifile,outdir, varargin)
% Prerequisites:
%   CanlabCore: 
%   yangchogosu_toolbox: https://github.com/didch1789/yanchogosu_toolbox.git
%   FSL
%
% inputs:
%       -inputdir: RAW nifti file.
%       -outdir: '/wanted/outcome/file/path'
%   varargins
%       -disdaq_num: number of imgs for disdaq (integer)
%       -TR: repetition time
%       -doMC: do motion correction(output filename contains ~moco. no following varargin) 
%       -duration: duration in seconds
%       -trial_dat: 
%       -nodisplay: Don't use fsleyes for final result.
%                   Display will show your beta map with the range within 95~99.9 percent.
%
% Primarily run in condition of DISDAQ and MCF (2020.07.16)
% e.g)
% inputdir = '/some/where/in/NAS/' or '/just/folder/that/contains/dicm/file'
% outputdir = 'some/where/in/your/local/drive'
% funclocal({inputdir}, outputdir, 'TR', 2, 'duration', 8, 'doMC')
%
if ~exist(outdir, 'dir'),mkdir(outdir);end

doMC = 0;
doJJ = 0;
displayFSL = true;

for i = 1:numel(varargin)
    switch varargin{i}
        case 'disdaq_num'
            disdaq_num = varargin{i+1};
        case 'TR'
            TR = varargin{i+1};
        case 'duration'
            dur = varargin{i+1};
        case 'doMC'
            doMC = 1;
        case 'doJJ' % special varagin only for LSD_Pain!!!
            doJJ = 1;
        case 'trial_dat'
            trialdat = varargin{i+1};
        case 'nodisplay'
            displayFSL = false;
    end
end

% 1.Disdaq
niifile = strrep(niifile, ' ' , '_');
[~, Fname, ext] = fileparts(niifile);
cmd1 = 'export FSLOUTPUTTYPE=NIFTI';

if exist('disdaq_num', 'var')
    disp('DO DISDAQ!')
    wholetr = numel(spm_vol(niifile));
    cmd2 = sprintf('fslroi %s %s %d %d', niifile, fullfile(outdir, [Fname, '_disdaq', ext]), ...
       disdaq_num , wholetr - disdaq_num);
    system([cmd1 ';' cmd2])
else
    disp('NO DISDAQ!')
end

% 2. BET on run1 disdaq
disp('BET on run1 disdaq!!!')
disdaqimg = sort_ycgosu(fullfile(outdir, '*disdaq*'));
disdaqimg = strrep(disdaqimg, ' ', '_');
cmd2 = sprintf('bet %s %s', disdaqimg{1}, fullfile(outdir, 'bet_with_disdaq.nii'));
system([cmd1,';',cmd2]);

% 4. FSL merge
% disp('Merging disdaq dats (only 1st img!!!)')
% disdaqs = sort_ycgosu([outdir, filesep, 'RUN*TR*', filesep, '*_disdaq.nii']);
% disdaqs = strrep(disdaqs, ' ', '_');
% cmdstr = ['fslmerge -tr ' fullfile(outdir, 'RUN_disdaq_concat.nii')];

% for i = 1:numel(disdaqs)
%     cmdstr = [cmdstr, blanks(1), disdaqs{i}];
% end
% 
% clipboard('copy', [cmd1, ';', cmdstr, blanks(1), num2str(TR)]);
% system(clipboard('paste'));

% 5. Motion correction!
if doMC
    disp('Motion correction-ing in disdaq imgs. 1st img is the reference volume!!!')
    cmd2 = ['mcflirt -in ', disdaqimg{1}, ...
        ' -o ',  fullfile(outdir, 'run_disdaq_moco.nii'),' -refvol 0 -plots'];
    system([cmd1, ';' cmd2]);
else
    disp('No motion correction!!!')
end 

% 6. GLM
if ~doJJ
    onsets = [];
    for i = 1:numel(trialdat)
        load(trialdat{i}) % data
        for j = 1:numel(data.dat.trial_dat)
            if exist('disdaq_num', 'var')
                onsets(j, i) = data.dat.trial_dat{j}.stim_starttime - data.dat.experiment_starttime;
            else
                onsets(j, i) = data.dat.trial_dat{j}.stim_starttime - data.dat.runscan_starttime;
            end
        end 
    end
else
    onsets = [0:25:334-25]';
end

Dm = [];
for i = 1:size(onsets, 2)
    reg_temp = onsets2fmridesign([onsets(:, i), repmat(dur, size(onsets, 1), 1)], ...
        TR, TR*(wholetr - disdaq_num), spm_hrf(TR));
    Dm(:, i) = reg_temp(:, 1);
end
Dm = Dm(:);

cd(outdir)
A = fmri_data('run_disdaq_moco.nii', 'bet_with_disdaq.nii');
mvmt = importdata('run_disdaq_moco.nii.par');
A.X = [Dm, mvmt, zscore(1:size(Dm, 1))']; % add linear detrend but didn't make a big difference!
out = regress(A, 'nodisplay');

betamap = out.b;
betamap.p(:, 2:end) = [];
betamap.ste(:, 2:end) = [];
betamap.threshold(:, 2:end) = [];
betamap.sig(:, 2:end) = [];
betamap.dat(:, 2:end) = [];
betamap.fullpath = fullfile(outdir, 'betamap.nii');
betamap.write;

disp('DONE!!!')
a1 = prctile(abs(betamap.dat), 99);
a2 = prctile(abs(betamap.dat), 99.9);
fprintf('Values are around %f ~ %f\n', a1, a2)

if displayFSL
    system(sprintf("fsleyes bet_with_disdaq.nii betamap.nii -dr %.2f %.2f -cm hot&", a1, a2))
end

end






