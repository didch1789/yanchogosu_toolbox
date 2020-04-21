%% Run on local computer
[rootdir, basedir, gitdir] = set_path_env('ncpu', 1);

%% Run on server
[rootdir, basedir, gitdir] = set_path_env('ncpu', 1, 'hostname', 'cnir9');

%% 1:4 5:8 9:12 13:16 17:20 21:24 25:28 29:32 33:36 37:40 41:44 45:48 49:52
sj_num = {51 52}; 
task_name = 'REST';
div_num = 10;
threshold = 0.05;
run_preproc_command = [];
for i = 1:numel(sj_num)
    run_preproc_command = [run_preproc_command 'matlab_orig -nodesktop -nosplash -nodisplay -r ' ...
        '"cd ' fullfile(basedir, 'projects/CAPS2/scripts') '; CAPS2_get_network_measures([' num2str(sj_num{i}) '], ''task_name'', ''' task_name ''', ''div_num'', ' num2str(div_num) ', ''threshold'', ' num2str(threshold) '); quit"' ...
        ' >> ~/' sprintf('sub-caps%.3d_network_log.txt', i) ' 2>&1 < /dev/null &'];
end
clipboard('copy', run_preproc_command);

%%
sj_num = 47:52;
run_preproc_command = [];
for i = 1:numel(sj_num)
    run_preproc_command = [run_preproc_command 'matlab_orig -nodesktop -nosplash -nodisplay -r ' ...
        '"cd ' fullfile(basedir, 'projects/CAPS2/scripts') '; CAPS2_preproc_network([' num2str(sj_num(i)) ']); quit"' ...
        ' >> ~/' sprintf('sub-caps%.3d_preproc_log.txt', sj_num(i)) ' 2>&1 < /dev/null &'];
end
clipboard('copy', run_preproc_command);

%% 1:4 5:8 9:12 13:16 17:20 21:24 25:28 29:32 33:36 37:40 41:44 45:48 49:52
sj_num = 52;
task_name = 'REST';
div_num = 1;
threshold = 0.10;
run_preproc_command = [];
for i = 1:numel(sj_num)
    run_preproc_command = [run_preproc_command 'nohup matlab_orig -nodesktop -nosplash -nodisplay -r ' ...
        '"cd ' fullfile(basedir, 'projects/CAPS2/scripts') '; CAPS2_get_network_gradients([' num2str(sj_num(i)) '], ''task_name'', ''' task_name ''', ''div_num'', ' num2str(div_num) ', ''threshold'', ' num2str(threshold) '); quit"' ...
        ' >> ~/' sprintf('sub-caps%.3d_task_%s_network_gradient_log.txt', sj_num(i), task_name) ' 2>&1 < /dev/null &'];
end
clipboard('copy', run_preproc_command);

%% 
sj_num = [1:52];
task_name = 'REST';
div_num = 10;
threshold = 0.05;
repeat_num = 100;
run_preproc_command = [];
for i = 1:numel(sj_num)
    run_preproc_command = [run_preproc_command 'nohup matlab_orig -nodesktop -nosplash -nodisplay -r ' ...
        '"cd ' fullfile(basedir, 'projects/CAPS2/scripts') '; CAPS2_get_module([' num2str(sj_num(i)) '], ''task_name'', ''' task_name ''', ''div_num'', ' num2str(div_num) ', ''threshold'', ' num2str(threshold) ', ''repeat_num'', ' num2str(repeat_num) '); quit"' ...
        ' >> ~/' sprintf('sub-caps%.3d_task_%s_network_module_log.txt', sj_num(i), task_name) ' 2>&1 < /dev/null &'];
end
clipboard('copy', run_preproc_command);

%%
div_idx = 9:10;
run_preproc_command = [];
for i = 1:numel(div_idx)
    run_preproc_command = [run_preproc_command 'nohup matlab_orig -nodesktop -nosplash -nodisplay -r ' ...
        '"cd ' fullfile(basedir, 'projects/CAPS2/scripts') '; CAPS2_community_contingency_multilayer_v2([' num2str(div_idx(i)) ']); quit"' ...
        ' >> ~/' sprintf('sub-caps%.3d_network_allegiance_log.txt', div_idx(i)) ' 2>&1 < /dev/null &'];
end
clipboard('copy', run_preproc_command);

system(run_preproc_command)

system('cd ~')



