% Basic directory setting 

md = 7;

% basedir = '/cocoanlab/habenula/';
% datdir = fullfile(basedir, 'projects/CAPS_project/data');
% pdatdir = fullfile(basedir, 'data/CAPS_data/conn_procdata');
% 
% Resource_dir = '/cocoanlab/habenula/Resources_server';
% addpath(genpath(Resource_dir));

basedir = '/Volumes/habenula/hbmnas/';
datdir = fullfile(basedir, 'projects/CAPS_project/data');
pdatdir = fullfile(basedir, 'data/CAPS_data/conn_procdata');

spider_dir = '/Resources/github/canlab/CanlabCore/CanlabCore/External/spider';
addpath(genpath(spider_dir));

% basedir = '/Users/jaejoonglee/hbmnas/habenula/';
% datdir = fullfile(basedir, 'projects/CAPS_project/data');
% pdatdir = fullfile(basedir, 'data/CAPS_data/conn_procdata');
% 
% spider_dir = '/Users/jaejoonglee/Documents/spider';
% addpath(genpath(spider_dir));

savedir = fullfile(datdir, ['conn_data/CAPS_dynamic_DCC_chunk_CPM_numc_model', strrep(num2str(md), ' ', ''), '.mat']);

%% Load data
load(fullfile(datdir, 'CAPS_dataset_NaN_170604.mat'));
load(fullfile(datdir, 'CAPS_painrating_170412.mat'));
load(fullfile(pdatdir, ['CAPS_DATA_whiten_total_model', strrep(num2str(md), ' ', ''), '_170922.mat']));

%% Indexing
good_subj = D.Subj_Level.data(:,6)==1;
nsubj = sum(good_subj==1);

conn_model = {
    'CAPS{subj_i}.roi_vals{1}.dat', ... % Buckner. 475 regions
    'CAPS{subj_i}.roi_vals{2}.dat', ... % Fan. 273 regions
    'CAPS{subj_i}.pexp_vals{1}.dat', ... % NPS_positive. 8 regions
    'CAPS{subj_i}.pexp_vals{2}.dat', ... % NPS_negative. 7 regions
    'CAPS{subj_i}.pexp_vals{3}.dat', ... % SIIPS1. 55 regions
    '[CAPS{subj_i}.pexp_vals{1}.dat, CAPS{subj_i}.pexp_vals{2}.dat, CAPS{subj_i}.pexp_vals{3}.dat]', ... % NPS_pos/NPS_neg/SIIPS1. total 70 regions
    'CAPS{subj_i}.roi_vals{3}.dat', ... % Fan + Brainstem. 279 regions
    'CAPS{subj_i}.roi_vals{4}.dat', ... % Glasser + Subcortex. 249 region
    };


TR = 0.46;

for j = 1:7
    if j == 1
        cidx{j} = 1 : round(45/TR*j - 20);
    else
        cidx{j} = round(45/TR*(j-1) - 20)+1 : round(45/TR*j - 20);
    end
end


temp_caps = capsaicin;
temp_rest = resting;

clear capsaicin resting;

capsaicin.intensity = cat(2,temp_caps.intensity{good_subj});
resting.intensity = cat(2,temp_rest.intensity{good_subj});
capsaicin.unpleasant = cat(2,temp_caps.unpleasant{good_subj});
resting.unpleasant = cat(2,temp_rest.unpleasant{good_subj});



%% Getting connectivity prediction

pred_model = {};

for model_i = md
    
    wh_fold = [];
    
    dat = fmri_data;
    dat.dat = [];
    
    for subj_i = 1:nsubj
        for chunk_i = 1:numel(cidx)
            wh_fold = [wh_fold; repmat(subj_i, 2, 1)];
            dat.dat = [dat.dat, mean(flat_cont_r{subj_i}{model_i}(:, cidx{chunk_i}), 2), ...
                mean(flat_caps_r{subj_i}{model_i}(:, cidx{chunk_i}), 2)];
            dat.Y = [dat.Y; resting.intensity(chunk_i, subj_i); capsaicin.intensity(chunk_i, subj_i)];
        end
    end
    
    numc = [0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001, 0.000005, 0.000001, ...
        0.0000005, 0.0000001, 0.00000005, 0.00000001, 0.000000005, 0.000000001, 0.0000000005, 0.0000000001];
    
    pred_model{model_i} = predict_CAPS_multiple(dat, 'algorithm_name', 'cv_cpm_numc', ...
        'polarity', 'both', 'threshnum', numc, 'polynorm', 1, 'spearman', 'nfolds', wh_fold, 'error_type', 'mse'); %0.000001
    
end



%% get prediction-outcome correlation
pred_out_corr = {};
pred_out_corr_mean = [];

for model_i = md
    for numc_i = 1:numel(numc)
        
        M = [];
        
        for subj_i = 1:nsubj
            pidx = numel(cidx)*2*(subj_i-1)+1 : numel(cidx)*2*(subj_i);
            pred_out_corr{model_i}{numc_i}{subj_i} = corr(pred_model{model_i}{numc_i}.stats.Y(pidx), pred_model{model_i}{numc_i}.stats.yfit(pidx));
            M = [M pred_out_corr{model_i}{numc_i}{subj_i}];
        end
        
        pred_out_corr_mean(model_i, numc_i) = mean(M);
        
    end
end

%plot(numc, pred_out_corr_mean);


%% save files
        
save(savedir, 'pred_*', 'numc')  
