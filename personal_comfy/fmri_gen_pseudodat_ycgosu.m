function out = fmri_gen_pseudodat_ycgosu(numtr, numregs, varargin)
% 
% pseudo random fmri_data with AR(1) random noise

out = zeros(numtr, numregs);
for i = 1:numel(varargin)
    switch varargin{i}
        case 'noisetype'
            noisetype = varargin{i+1};
    end
end            

for i = 1:numregs
    x = randi([1, 100], 1, 2);
    y = randi([1 3], 1, 1);
    z = randi([1 5], 1, 1);
    
    if strcmp(noisetype, 'ar')
        Mdl = arima('Constant', 0, 'AR', {0.5}, 'Variance', z);
        noise = simulate(Mdl, numtr);
    elseif strcmp(noisetype, 'gaussian')
        noise = randn(numtr, 1);
    end
    
    xi = y*sin(linspace(x(1), x(2), numtr))' + z*noise;
    out(:, i) = xi;
end

