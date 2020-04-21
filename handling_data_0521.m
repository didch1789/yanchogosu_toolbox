%% variance value

data = fmri_data;
data.dat = [];

var_dcc = [];
sub_idx = find(dat.Both_idx);

for i = 1:numel(sub_idx)
    var_dcc(:, i) = var(dat.Dcc{sub_idx(i)},0,2);
end

%% mean value 

mean_dcc = [];

for i = 1:numel(sub_idx);
    mean_dcc(:, i) = mean(dat.Dcc{sub_idx(i)}, 2);
end

%% leave-one-out prediction using variance

data2 = data; 
data2.dat = scale(var_dcc);

data2.Y = scale(dat.Behavioral.RRS_BroodRum(sub_idx));

[cverr, stats, optout] = predict(data2, 'algorithm_name', 'cv_pcr', 'nfolds', numel(data.Y), 'error_type', 'mse', 'numcomponents', 4);

stats_var = stats;

% r =  -0.3087 
% using boxcox also yields -0.2687

%% leave-one-out prediction using mean

data2.dat = scale(mean_dcc);

data2.Y = scale(dat.Behavioral.RRS_BroodRum(sub_idx));

[cverr, stats, optout] = predict(data2, 'algorithm_name', 'cv_pcr', 'nfolds', numel(data.Y), 'error_type', 'mse', 'numcomponents', 4);

stats_mean = stats;

% r = 0.0389

%% try other dependent variables(RRS_Rum)

data3 = fmri_data;

data3.Y = scale(dat.Behavioral.RRS_Rumination(sub_idx));

data3.dat = scale(var_dcc);
%data3.dat = scale(mean_dcc)

[cverr, stats, optout] = predict(data3, 'algorithm_name', 'cv_pcr', 'nfolds', 1:numel(data.Y), 'error_type', 'mse');

stats_rum = stats;

% for mean : r = -0.1072
% for var : r = -0.2922

%% try other dependent variables(RRS_Depression)

data3.Y = scale(dat.Behavioral.RRS_Depression(sub_idx));

data3.dat = scale(mean_dcc);
%data3.dat = scale(var_dcc);

[cverr, stats, optout] = predict(data3, 'algorithm_name', 'cv_pcr', 'nfolds', 1:numel(data.Y), 'error_type', 'mse');

stats_dep = stats;

% r = 0.0061

%% try other algorithm(lasso_pcr)
data3 = fmri_data;

data3.Y = scale(dat.Behavioral.RRS_Depression(sub_idx));
data3.dat = scale(var_dcc);

[cverr, stats, optout] = predict(data3, 'algorithm_name', 'cv_lassopcr', 'lasso_num', 60, 'nfolds', 1:numel(data3.Y), 'error_type', 'mse', ...
    'EstimatedParams');

% r : -0.2503



%% pcr_numc (jj)

addpath('/Volumes/wissen/cocoanlab Dropbox/projects/bmrk5/script');
data4 = fmri_data;
data4.Y = scale(dat.Behavioral.RRS_BroodRum(dat.Both_idx));
data4.dat = mean_dcc;

numc = 1:80;
% number of principal components to run.
pred_model = predict_CAPS_multiple(data4, 'algorithm_name', 'cv_pcr_numc', 'numcomponents', numc, 'nfolds', 1:numel(data4.Y), 'error_type', 'mse');

%% cpm_numc

numc = [0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001];
pred_model = predict_CAPS_multiple(data4, 'algorithm_name', 'cv_cpm_numc', ...
    'polarity', 'both', 'threshnum', numc, 'polynorm', 1, 'spearman', 'nfolds', 1:numel(data4.Y), 'error_type', 'mse');

%% relation between dcc_mean and dcc_variance
corrs = [];

for i = 1:numel(sub_idx)
    corrs(:, i) = corr(mean_dcc(:, i), var_dcc(:, i));
end

% cuz of skewd data.Y ?