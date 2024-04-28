function make_batch_script_ycgosu(funcdir, funcname, logdir, batchidx)
% funcdir: should be in ''' ''''. 
%   e.g) '''/your/path/'''
% funcname: function name without .m in 'funcname'
% logdir: directory of log directory in 'log directory'
% numbatch: # of batch


if ~exist(logdir, 'dir')
    mkdir(logdir)
end

clearvars s
s = '';
randstr = randseq_generator(5);
for i=batchidx
 x = ['nohup matlab -singleCompThread -batch "cd(', funcdir, ');', funcname, '(', num2str(i),  ')"'];
 y = [' -logfile ', fullfile(logdir, sprintf('log%02d-%s.txt', i, randstr)), '&', ' '];
 s = [s, x, y];
end

clipboard('copy', s)

end