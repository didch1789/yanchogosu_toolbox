%% mean value 

mean_dcc = [];

for i = 1:numel(sub_idx);
    mean_dcc(:, i) = mean(dat.Dcc{sub_idx(i)}, 2);
end