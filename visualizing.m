%% load data
basedir = '/Users/jungwookim/Dropbox/Labstuff/BMRK5/BMRK5_matlab/';
uiopen(fullfile(basedir,'bmrk5_compiled_measures_161012_UPDATED.csv'), 1);
    %should be a better way...
compiled = bmrk5compiledmeasures161012UPDATED;

variables_of_interest = {'Subject_ID', 'Gender', 'Ethnicity', 'Age',...
    'RSFC_Pos', 'RSFC_Neg', 'RSFC_Centrality', 'RSFC_Social',...
    'RSFC_Imagery', 'RSFC_Present', 'RSFC_Past', 'RSFC_Future', 'BDI_Sum', 'RRS_Brood', 'RRS_Depression', ...
    'RRS_Rumination', 'RRS_Sum', 'RRS_BroodRum', 'Avg_Imagery', 'Avg_Pos', 'Avg_Neg', 'Avg_Control', 'Avg_Centrality', ...
    'Avg_Social', 'Avg_Duration', 'Avg_Durrev', 'FA_Val_NegMean', 'FA_Val_PosMean', 'FA_Val_NeutMean', 'FA_Val_TotMean', ...
    'FA_Rep_NegSum', 'FA_Rep_NeutSum', 'FA_Rep_PosSum', 'FA_LSA_NegMean', 'FA_LSA_NeuMean',...
    'FA_LSA_PosMean', 'FA_LSA_TotMean'};

compiled_table = compiled(:, variables_of_interest);

compiled_table1 = compiled_table(sum(ismissing(compiled_table, NaN), 2) < 8, :);
    % subject = 88
    % compiled_table2 = compiled_table(sum(ismissing(compiled_table, NaN), 2) < 9, :);
        % subject = 92
        % 891, 983, 1023, 1024, 1114
compiled_cell = table2cell(compiled_table1);
compiled_measures = table2struct(compiled_table1, 'ToScalar', true);

%%
load(which('cluster_Fan_Net_r280.mat'));
cluster_idx = cluster_Fan_Net.dat(:,3);
% see description

cols = [ 0 0.4470 0.7410
 0.8500 0.3250 0.0980
 0.9290 0.6940 0.1250
 0.4940 0.1840 0.5560
 0.4660 0.6740 0.1880
 0.3010 0.7450 0.9330
 0.6350 0.0780 0.1840
 0.9843 0.6039 0.6000
 0.7922 0.6980 0.8392];
% color value

%% 
load('dat_obj_bmrk5_S998_OC1465_preproc_Brainnetome_29-Mar-2019.mat', 'dcc_val');
dcc_val_998 = dcc_val;

var_1192=var(dcc_val_998');
reform_1192_var = reformat_r_new(var_1192, 'reconstruct');


vis_corr(reform_1192_var, 'colorbar', 'nolines', 'group', cluster_idx, 'group_linewidth', 5, 'group_linecolor', 'k', ...
 'group_tick', 'group_tickstyle', 'edge', 'group_tickwidth', 2, 'group_ticklength', 5, 'group_tickoffset', 3, 'smooth',...
 'clim', [0.04 0.06], 'triangle', 'triangle_color', [0 0.4392 0.7529], 'triangle_width', 7);
set(gcf, 'position', [1000  745  682  593]);


out_rest = vis_corr(reform_998_var, 'colorbar', 'nolines', 'group', cluster_idx, 'group_mean', 'group_color', cols,...
    'triangle', 'triangle_color', [0 0.4392 0.7529], 'triangle_width', 5);
set(gcf, 'position', [1000  745  682  593]);
