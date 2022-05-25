function aIdx = generate_randAligns_ycgosu(mat1, mat2, varargin)

%% This function is to generate random alignment index. 
% Basic principle is that,
% if there is a random data that share same covariance structure with
% the real data, will that random data show orthogonality?
n_dim = 10;
for i = 1:numel(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'n_dim'
                n_dim = varargin{i+1};
        end
    end
end
all_epoch_data = [mat1;mat2];
n_time = size(all_epoch_data, 1);
n_boot = 1000;

aIdx = struct;
for iB = 1:n_boot
    rng(iB + n_dim);
    perm_epoch_data = all_epoch_data(randperm(n_time), :);
    % permutation do not change covariance.
    perm_mat1 = perm_epoch_data(1:size(mat1, 1), :);
    perm_mat2 = perm_epoch_data(size(mat1, 1)+1:end, :);
    [mat1PCs, ~, mat1Var] = pca(perm_mat1, 'NumComponents', n_dim);
    [mat2PCs, ~, mat2Var] = pca(perm_mat2, 'NumComponents', n_dim);

%     expl1on2 = trace(mat2PCs' * cov(perm_mat1) * mat2PCs) / trace(mat1PCs' * cov(perm_mat1) * mat1PCs);
%     expl2on1 = trace(mat1PCs' * cov(perm_mat2) * mat1PCs) / trace(mat2PCs' * cov(perm_mat2) * mat2PCs);
    expl1on2 = sum(var((perm_mat1 - mean(perm_mat1)) * mat2PCs)) ./ sum(mat1Var(1:n_dim));
    expl2on1 = sum(var((perm_mat2 - mean(perm_mat2)) * mat1PCs)) ./ sum(mat2Var(1:n_dim));
    % cov takes a while if #feauture is large.
    
    aIdx.mean(iB, :)     = (expl1on2 + expl2on1)/2;
    aIdx.expl1on2(iB, :) = expl1on2;
    aIdx.expl2on1(iB, :) = expl2on1;
end
    
	
end