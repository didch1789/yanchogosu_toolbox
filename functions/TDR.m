function out = TDR(Xin, Yin, varargin)
% Input:
%   Xin: regressors encoding task variables(t.v.). 
%        Size of #trial X #t.v.
%        (e.g., 1 task variable (e.g., amount of reward) can be encoded as
%         series of ones and minus ones. 1 for high reward and -1 for low reward.)
%        NOTES1. scale of each task variable should be (ideally) same, to avoid
%        variable specific over(or under)estimation of regression coefficient.
%        NOTES2. This matrix also determines number of condition.
%                If the matrix has the following values,
%                   [1 0;
%                    0 1;
%                    0 0;
%                    1 0;
%                    0 1;
%                    0 0];
%                It is regarded as 6 trials with 3 conditions and 2 task
%                variables.
%   Yin: neural firing rates of trials.
%        Size of #time X #trial X #neuron
%   varargin:
%       ('norm_y', boolean) : zscoring of each neuron. (defaults: true)
%                  (mean and std are computed after combining
%                   time and trials.)
%       ('intcpt', boolean)    : Whether to add intercept in the GLM process. (defaults: true)
%       ('interact', boolean)  : Whether to add interaction term. (defaults: false)
%                                This will add interaction term of t.v.s.
%       ('binsize', integer)   : Definition of bin size. 
%                                If not specified, it won't do any binning.
%       ('smoothsize', integer): Definition of smoothing window. (defaults: 30)
%       ('numcomp', integer)  : Definition of number of PCs. 
%                               (defaults: 90% explaining PCs).
%
% Output:
%   out.regression_weights : #neuron X #t.v. X #time
%                            (encoding weights in time.)
%   out.beta_v_orth : #neuron X #t.v.
%                    (t.v. specific axis)
%   out.p_vc: #time X #conds X #t.v.  
%             (average population rates projected. 
%   #time can differ in "Ouput" variable from "Input" variable
%   depending on "binsize".
%   
%  (Example)
%  rng('default')
%  egXin = randi([0 1], 17, 2);
%  egYin = randn(100, 17, 20); 
%  out = TDR(egXin, egYin, 'binsize', 10, 'smoothsize', 20, 'numcomp', 5);
%  out = 
%     regression_weights: [20×2×10 double]
%            beta_v_orth: [20×2 double]
%                   p_vc: [10×4×2 double]
%
%
%   2022.03.06. Jungwoo.

norm_y = true;
add_intcpt = true; 
add_interact = false;
smthwindow = 30;

for i = 1:numel(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'norm_y', 'normy'}
                % zscoring neural firing rates.
                norm_y = varargin{i+1};
            case {'intcpt'}
                add_intcpt = varargin{i+1};
            case {'interact'}
                add_interact = varargin{i+1};
            case {'binsize'}
                binsize = varargin{i+1};
            case {'smoothsize'}
                smthwindow = varargin{i+1};
            case {'numcomp'}
                numcomp = varargin{i+1};
        end
   end
end

origY = Yin;
origX = Xin;
[n_time, n_trial, n_neuron] = size(origY);
[n_trial2, n_tv] = size(origX);
if n_trial2 ~= n_trial
    error('Number of trials in X and Y are not equal.')
end

% Specifying Y.
Y = reshape(Yin, [], n_neuron);
if exist('binsize', 'var')
    if mod(n_time, binsize) == 0
        Y = bin_raster(Y', binsize)';
    else
        error('Choose appropriate binsize!')
    end
end
if norm_y
    Y = zscore(Y);
    Y = reshape(Y, [], n_trial, n_neuron);
end

% Specifying X.
if add_interact
    inter_var = nchoosek(1:n_tv, 2);
    inter_tvs = NaN(size(Xin, 1), size(inter_var, 1)); % prealloc.
    for ii = 1:size(inter_var, 1)
        inter_tvs(:, ii) = Xin(:, inter_var(ii, 1)) .* Xin(:, inter_var(ii, 2));
    end
    Xin = [Xin inter_tvs];
end
if add_intcpt
    Xin = [Xin ones(size(Xin, 1), 1)];
end
X = Xin;

%% Regress Y with respect to X.    
glmweights = NaN(n_neuron, n_tv, size(Y, 1)); % prealloc. 
for ni = 1:n_neuron
    trXti = Y(:, :, ni)';
    glmweight = pinv(X) * trXti;
    glmweights(ni, :, :) = glmweight(1:n_tv, :);
end
out.regression_weights = glmweights(:, 1:n_tv, :); 


%% PCA (6.4 - 6.7)
[~, ~, condIDs] = unique(origX, 'rows');
n_cond = numel(unique(condIDs));

Yconds = [];
for ci = 1:n_cond
    Yavgsmth = smoothdata(squeeze(mean(Y(:, condIDs == ci, :), 2)), 1, 'gaussian', smthwindow);
    Yconds = cat(1, Yconds, Yavgsmth);
end

[PCcoeff, ~, ~, ~, expls] = pca(Yconds);
if ~exist('numcomp', 'var')
    numcomp = find(cumsum(expls) >=90, 1);
end

PCcoeff = PCcoeff(:, 1:numcomp);
D = PCcoeff * PCcoeff';

B_max = NaN(n_neuron, n_tv);
for si = 1:n_tv
    beta_vt = squeeze(glmweights(:, si, :));
    beta_vt_PC = D*beta_vt;
    [~, maxidx] = max(sum(beta_vt_PC .^ 2, 1));
    B_max(:, si) = beta_vt_PC(:, maxidx).^2;
end

[Q, ~] = qr(B_max);
%  Q is calculated based on first column of B_max.
%  corr(Q(:, 1), B_max(:, 1)) results in  '-1 '
beta_v_orth = Q(:, 1:n_tv);
% beta_v_orth determines t.v. specific axis.
p_vc = Yconds * beta_v_orth;
p_vc = reshape(p_vc, size(Y, 1), [], n_tv);
% p_vc now has size of #time X #conds X #t.v. 
% p_vc(:, :, 1) which has size of #time X #conds indicates
% average firing rates projected on to axis which contains 
% variable specifi information.

out.beta_v_orth = beta_v_orth;
out.p_vc = p_vc;

%% sub-function
function binned = bin_raster(rasterdat, binsize)
% Inputs
%   rasterdat : raw spike data in 1ms. (neurons X timeseries)
%   binsize   : time bin size in ms (e.g. 10) 
%
% Outputs
%   binned: binned (summed wihthin bin) raster 
% e.g. 
%   size(rasterdat) = 100 200 (100 neurons' 200 ms data)
%   binsize = 10;
%   size(binned)    = 100 20
%
%   size(rasterdat) = 100 205 (100 neurons' 200 ms data)
%   binsize = 10
%   size(binned)    = 100 20 (first 5 bin will contain 11 binsize)
neuronnum   = size(rasterdat, 1); 
timesize    = size(rasterdat, 2);
Quotient    = fix(timesize / binsize);
R           = mod(timesize, binsize);
binned      = NaN(neuronnum, Quotient);

l = 1;
for bb = 1:Quotient
    if R ~= 0
        k = binsize + 1;
        R = R - 1;
    else
        k = binsize;
    end
    
    binned(:, bb) = sum(rasterdat(:, l:l+k-1), 2);
    l = l + k;

end

end

end