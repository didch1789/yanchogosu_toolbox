function out = fmri_gen_pseudodat_ycgosu(numtr, numregs)
% pseudo random fmri_data with AR(1) random noise

out = randn(numtr, numregs);

for i = 1:numregs
    x = randi([0, 100], 1, 2);
    y = randi([-4 4], 1, 1);
    z = randi([1 3], 1, 1);
    Mdl = arima('Constant', 0, 'AR', {0.5}, 'Variance', z);
    noise = simulate(Mdl, numtr);
    xi = y*sin(linspace(x(1), x(2), numtr))' + z*noise;
    out(:, i) = xi;
end

