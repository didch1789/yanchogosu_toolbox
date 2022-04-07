function out = calc_trial_num(condIds)
n_sess = numel(condIds);
for i = 1:n_sess
    sess_cond = condIds{i};
    nanidx = isnan(condIds{i});
    int_condId = sess_cond(~nanidx);
    [cond_digits, ~, cond_Ids] = unique(int_condId);

    for j = 1:numel(cond_digits)
        out(j, i) = sum(cond_Ids == j);
    end
end

end