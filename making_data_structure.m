%% load data
basedir = '/Volumes/wissen/cocoanlab Dropbox/projects/bmrk5/sync/data/';
M = readtable(fullfile(basedir,'Behavioral_data.csv'));
    % when csv file contains both numeric and text value, readtable seems
    % like a good option
M.Properties.VariableNames{1} = 'Subject_ID';
    % idk why but 'Subject_ID' variable name gets changed when it is imported...
    
variables_of_interest = {'Subject_ID', 'Gender', 'Ethnicity', 'Age',...
    'RSFC_Pos', 'RSFC_Neg', 'RSFC_Centrality', 'RSFC_Social',...
    'RSFC_Imagery', 'RSFC_Present', 'RSFC_Past', 'RSFC_Future', 'BDI_Sum', 'RRS_Brood', 'RRS_Depression', ...
    'RRS_Rumination', 'RRS_Sum', 'RRS_BroodRum', 'Avg_Imagery', 'Avg_Pos', 'Avg_Neg', 'Avg_Control', 'Avg_Centrality', ...
    'Avg_Social', 'Avg_Duration', 'Avg_Durrev', 'FA_Val_NegMean', 'FA_Val_PosMean', 'FA_Val_NeutMean', 'FA_Val_TotMean', ...
    'FA_Rep_NegSum', 'FA_Rep_NeutSum', 'FA_Rep_PosSum', 'FA_LSA_NegMean', 'FA_LSA_NeuMean',...
    'FA_LSA_PosMean', 'FA_LSA_TotMean'};

data_temp = M(:, variables_of_interest);
    % Think it's easier to arrange data when it's in table rather than in structure... 
    
%% Field1 : Subject_ID (cell array)

Subject_ID_temp = table2cell(data_temp(:, 1));
all_subs = table2array(data_temp(:, 1)) ;

%% Field2 : Gender (cell array)

Gender_temp = table2cell(data_temp(:, 2));

%% FIeld3 : Ethnicity (cell array)

Ethnicity_temp = table2cell(data_temp(:, 3));

%% Field4 : Age (cell array)

Age_temp = table2cell(data_temp(:, 4));

%% Field5 : Behavioral ( numeric array )

Behavioral_temp_temp = table2array(data_temp(:, 5:end));

% making Behavioral_temp data to class 'double' and putting them into array
% and changing '.' to NaN
Behavioral_temp2 = [];
for i = 1:size(Behavioral_temp_temp, 1)
    for j = 1:size(Behavioral_temp_temp, 2)
        if class(Behavioral_temp_temp{i , j}) == 'char'
            try
                Behavioral_temp2(i, j) = str2num(Behavioral_temp_temp{i , j});
            catch
                Behavioral_temp2(i, j) = NaN;
            end
        end
    end
end

% if using data in structure
Behavioral_temp_table = array2table(Behavioral_temp2, 'VariableNames', variables_of_interest(5:end));
Behavioral_temp = table2struct(Behavioral_temp_table, 'ToScalar', true);



%% Field6: DCC(cell array) / Brain_idx_temp is defined at Field8 
load('/Volumes/wissen/cocoanlab Dropbox/projects/bmrk5/sync/data/dcc_vals.mat');
dcc_names_temp = who('-file', '/Volumes/wissen/cocoanlab Dropbox/projects/bmrk5/sync/data/dcc_vals.mat');
temp_list = setdiff(all_subs, brain_null_idx);

% arranging filenames in ascending order
dcc_names_temp = sort(dcc_names_temp);
for i = 1:numel(dcc_names_temp)
    if length(dcc_names_temp{i}) ~= length(dcc_names_temp{i+1})
        break
    end
end

if length(dcc_names_temp{i}) > length(dcc_names_temp{i+1})
    filenames = cat(1, dcc_names_temp(i+1:end), dcc_names_temp(1:i));
else 
    filenames = cat(1, dcc_names_temp(1:i), dcc_names_temp(i+1:end));
end


% making NaN cell array
Dcc_temp = cell(110, 1);
for i = 1:length(Dcc_temp)
    if isempty(Dcc_temp{i})
        Dcc_temp{i} = NaN;
    end
end

% putting data into NaN cell array
j = 1;
for i = 1:numel(Dcc_temp)
    if Brain_idx_temp(i) == true
        eval(['Dcc_temp{' num2str(i) '} = ' filenames{j} ';']);
        j = j + 1;
    end
end
% now Dcc_temp has data in it.

% index checking...

% k = 0;
% for i = 1:numel(Dcc_temp)
%     if isnan(Dcc_temp{i})
%         k = k + 1;
%     end
% end
% k = 17


%% Field7 : Behav_idx ( logical array _ exist = 1 / ~exist = 0)
% missing some Avg_Duration, Avg_Durrev, BDI_SUM!!!!!
behav_null_idx = [876, 891, 967, 969, 976, 983, 999, 1003, 1019, 1020, 1023, 1024, ...
    1041, 1090, 1114, 1118, 1146, 1169, 1194, 1211, 1216, 1221];

m = [];
for i = 1:numel(all_subs)
    if ismember(str2double(all_subs{i}), behav_null_idx)
        m(i) = false;
    else
        m(i) = true;
    end
end

Behav_idx_temp = logical(m');

% sum(Behav_idx_temp) = 88
% among 88 [916, 923, 979, 1047]  brain data null

%% Field8: Brain_idx

m = [];
brain_null_idx = sort([916, 923, 967, 979, 983, 999, 1003, 1019, 1020, 1023, 1041, 1047, 1090, 1114, 1146, 1211, 1221]);

for i = 1:numel(all_subs)
    if ismember(all_subs(i), brain_null_idx)
        m(i) = false;
    else
        m(i) = true;
    end
end

Brain_idx_temp = logical(m');

% among 93 [1024 1118 1169 1194 1216 876 891 969 976] behav data null

%% Field9: Both_idx

m = [];
both_null_idx = [876, 891, 967, 969, 976, 983, 999, 1003, 1019, 1020, 1023, 1024, ...
    1041, 1090, 1114, 1118, 1146, 1169, 1194, 1211, 1216, 1221, 916, 923, 979, 1047];

for i = 1:numel(all_subs)
    if ismember(all_subs(i), both_null_idx)
        m(i) = false;
    else
        m(i) = true;
    end
end

Both_idx_temp = logical(m');

%% Make a structure!
fields_temp = {'Subject_ID_temp', 'Gender_temp', 'Ethnicity_temp', 'Age_temp', ...
    'Behavioral_temp', 'Dcc_temp', 'Behav_idx_temp', 'Brain_idx_temp', 'Both_idx_temp'};

for i = 1:numel(fields_temp)
    field_name_temp = split(fields_temp{i}, 'temp');
    field_name_temp2 = field_name_temp{1};
    field_name = field_name_temp2(1:end-1);
    eval(['dat.' field_name '=' fields_temp{i} ';']);
end
% structure 'dat' created
    
%% saving data
save('dat.mat','dat','-v7.3');
% for data larger than 2GB, additional arg required.
