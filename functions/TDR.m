function [out, Yout] = TDR(Xin, Yin, varargin)
% mainly based on supplementary information of Mante et al., Nature, 2013.
% Input:
%   Xin: regressors encoding task variables(t.v.). 
%        Size of #sessions x 1 cell array. Each session (cell) conatains #trial X
%        #t.v. array.
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
%        Size of #sessions x 1 cell array. Each session (cell) contains
%        #time X #trial X #neuron array.
%   varargin:
%       ('norm_y' or 'normy', boolean) : zscoring of each neuron. (defaults: true)
%                                       (mean and std are computed after combining
%                                        time and trials.)
%       ('intcpt', boolean)    : Whether to add intercept in the GLM process. (defaults: true)
%       ('interact', boolean)  : Whether to add interaction term. (defaults: false)
%                                This will add interaction term of t.v.s.
%       ('binsize', integer)   : Definition of bin size. 
%                                If not specified, it won't do any binning.
%       ('smoothsize', integer): Definition of smoothing window. (defaults: 30)
%       ('numcomp', integer)   : Definition of number of PCs. 
%                               (defaults: 90% explaining PCs).
%       ('mv_mean', integer)   : sliding window average firing rates in a specified time window
%                               If not specified, it won't do any move mean.

%  NOTES. One weird thing is that, "movmean" is used to calculate GLM
%  weights, while for pseudo-population, it uses "smoothing". But I think
%  there are no big differences between the two.
%  NOTES. Not sure whether binning is used... 
%
% Output:
%   out.regression_weights : #sess x 1 cell, each cell has size of #neuron X #t.v. X #time
%                            (encoding weights in time. This is calculated per session)
%                            
%   out.beta_v_orth : #neuron X #t.v.
%                    (t.v. specific axis)
%                     #neuron is number of neurons in all sessions.
%                     This is calculated by combining
%                     out.regression_weights (see line 172)
%   out.p_vc: #time X #conds X #t.v.  
%             (average population rates projected. 
%             #time can differ in "Ouput" variable from "Input" variable
%             depending on "binsize".
%   
%  (Example)
%  egXin = {randi([0 1], 17, 2)};
%  egYin = {randn(100, 17, 20)}; 
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
            case {'mv_mean', 'mvmean'}
                mv_window = varargin{i+1};
        end
   end
end

origY = Yin;
origX = Xin;

if numel(origX) ~= numel(origY)
    error('Number of sessions in X and Y are not equal.')
else
    n_sess = numel(origX);
    Yall = cell(n_sess, 1);
    Xall = cell(n_sess, 1);
end

for i = 1:n_sess 
    [n_time, n_trial, n_neuron] = size(origY{i});
    [n_trial2, n_tv] = size(origX{i});
    if n_trial2 ~= n_trial
        error('Number of trials in X and Y are not equal.')
    end

    % Specifying Y.
    Y = reshape(origY{i}, [], n_neuron);
    if exist('binsize', 'var')
        if mod(n_time, binsize) == 0
            Y = bin_raster(Y', binsize)';
        else
            error('Choose appropriate binsize!')
        end
    end
    Yout.Bin = reshape(Y, [], n_trial, n_neuron);
    
    if exist('mv_window', 'var')
        Yout.MvwBin = movmean(Yout.Bin, mv_window, 1);
%         Yout.MvwBin = smoothdata(Yout.Bin, 1, 'gaussian', mv_window);
    end
    
    if norm_y
        if contains('MvwBin', fields(Yout))
            Yout.ZscMvwBin = zscore(Yout.MvwBin, 0, 1:2);    
        else
            Yout.ZscBin = zscore(Yout.Bin, 0, 1:2);  
        end
    end
    Yall{i} = Yout;

    % Specifying X.
    X = origX{i};
    if add_interact
        inter_var = nchoosek(1:n_tv, 2);
        inter_tvs = NaN(size(X, 1), size(inter_var, 1)); % prealloc.
        for ii = 1:size(inter_var, 1)
            inter_tvs(:, ii) = X(:, inter_var(ii, 1)) .* X(:, inter_var(ii, 2));
        end
        X = [X inter_tvs];
    end
    if add_intcpt
        X = [X ones(size(X, 1), 1)];
    end
    Xall{i} = X;
end


%% Regress Y with respect to X.    
for i = 1:n_sess
    ymetrics = fields(Yall{i});
    Yint = Yall{i}.(ymetrics{end});
    n_neuron = size(Yint, 3);
    n_time   = size(Yint, 1);
    glmweights = NaN(n_neuron, n_tv, n_time); % prealloc. 
    for ni = 1:n_neuron
        trXti = Yint(:, :, ni)';
        glmweight = pinv(Xall{i}) * trXti;
        glmweights(ni, :, :) = glmweight(1:n_tv, :);
    end
    out.regression_weights{i} = glmweights(:, 1:n_tv, :); 
end
glmweights_all = cat(1, out.regression_weights{:});


%% PCA (6.4 - 6.7)
Ycondcell = {};
for i = 1:n_sess
    Yint = Yall{i}.Bin;
    [~, ~, condIDs] = unique(Xall{i}, 'rows');
    n_cond = numel(unique(condIDs));
    for ci = 1:n_cond
        Yavgsmth = smoothdata(squeeze(mean(Yint(:, condIDs == ci, :), 2)), 1, ...
            'gaussian', smthwindow);
        Ycondcell{ci, i} = Yavgsmth;
    end
end
Yout.SmtAvgBin = cell2mat(Ycondcell);
if norm_y
    Yout.ZscSmtAvgBin = zscore(Yout.SmtAvgBin);
end
n_neuron_all = size(Yout.SmtAvgBin, 2);
ymetrics = fields(Yout);
Ycondmat = Yout.(ymetrics{end});

[PCcoeff, ~, ~, ~, expls] = pca(Ycondmat);
if ~exist('numcomp', 'var')
    numcomp = find(cumsum(expls) >=90, 1);
end

PCcoeff = PCcoeff(:, 1:numcomp);
D = PCcoeff * PCcoeff';

B_max = NaN(n_neuron_all, n_tv);
for si = 1:n_tv
    beta_vt = squeeze(glmweights_all(:, si, :));
    beta_vt_PC = D*beta_vt;
    [~, maxidx] = max(sum(beta_vt_PC .^ 2, 1));
    B_max(:, si) = beta_vt_PC(:, maxidx); % stability of regression coefficients (rank order of beta_vt_PC.)
end

[Q, ~] = qr(B_max);
%  Q is calculated based on first column of B_max.
%  corr(Q(:, 1), B_max(:, 1)) results in  '-1 '
beta_v_orth = Q(:, 1:n_tv);
% beta_v_orth determines t.v. specific axis.
p_vc = Ycondmat * beta_v_orth;
p_vc = reshape(p_vc, [], n_cond, n_tv);
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