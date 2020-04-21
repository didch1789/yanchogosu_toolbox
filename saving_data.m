data_list = dir('dat_*.mat');

%% saving dcc as subject number

for i = 1:length(data_list)
    name_split = split(data_list(i).name, '_');
    subj_num = strcat(name_split{4}, '_', name_split{5});
    load(data_list(i).name, 'dcc_val');
    eval(['dcc_val_' subj_num ' = dcc_val;']);
    fprintf('%d th subject_finished, total: 93 \n', i)
end

save('dcc_vals.mat', 'dcc_val_*');


% broodrum upper 25%  {980,  877, 1089, 1125,  913,  979, 1217, 1109, 1219, 1113,  963,
%            1198,  947, 1040, 1077, 1086, 1161, 1078,  946, 1120, 1071, 1000}

% broodrum_25 = [980,  877, 1089, 1125,  913,  979, 1217, 1109, 1219, 1113,  963, 1198,...
%     947, 1040, 1077, 1086, 1161, 1078,  946, 1120, 1071, 1000];

%% saving_dcc_mean
for i = 1:length(data_list)
    name_split = split(data_list(i).name, '_');
    s_num = name_split{4};
    num = s_num(2:end);
    load(data_list(i).name, 'dcc_val')
    eval(['dcc_mean_' num '=mean(dcc_val,2);']);
end

save('dcc_mean_of_all.mat', '*mean*')

%% saving_dcc_variance
for i = 1:length(data_list)
    name_split = split(data_list(i).name, '_');
    s_num = name_split{4};
    num = s_num(2:end);
    load(data_list(i).name, 'dcc_val');
    dcc_var = var(dcc_val');
    eval(['dcc_var_' num '= dcc_var;']);
end

save('dcc_var_of_all.mat', '*var*')
    
    
    