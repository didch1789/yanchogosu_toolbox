function [stats, Xscores, Yscores, Xweights, Yweights] = predict_PLS_ycgosu(obj, varargin)

% Partial Least Squares cross-validation. (based on NIPALS algorithm or SIMPLS algorithm)
% inputs:
%   obj : fmri_data object; 
%       obj.dat: p x n matrix (p: features(voxels), n: oberservations)
%       obj.Y: n x 1 matrix
%   varargin : 
%       {'numcomponents', k} or {'numcomps', k} k: number of components want to
%                                              extract from PLS algorithm.
%       {'nfolds', q} q: number of cross validation. (default: leave-one-out)
% outputs:
%   stats:
%       beta : beta includes intercept at first row
%       weights
%       yfit
%       pred_outcome_r and p

% for SIMPLS no varargin will take maximum features it can.
% for NIPALS no varargin will automatically estimate the features.

Y = obj.Y;
X = obj.dat';


nY = size(Y, 1);
Xrow = size(X, 1);
Xcol = size(X, 2);

if Xrow ~= nY
    error('Size mismatch!!!')
end

ncomp = min(Xrow-1, Xcol); % default
nfolds = Xrow; % leave-one-out cv as defaults
cvpart = stratified_holdout_set(Y(:, 1), nfolds); % leave-one-out cv as default

if numel(varargin) ~= 0
    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}
                % functional commands
                case {'numcomponents', 'numcomps'}
                    ncomp = varargin{i+1};
                case {'nfolds'}
                    nfolds = varargin{i+1};
                    cvpart = stratified_holdout_set(Y, nfolds);
%                 case {'NIPALS', 'nipals'}
%                     nipals = 1;simpls = 0;
%                case {'SIMPLS', 'simpls'}
%                    simpls = 1;nipals = 0;
            end
        end
    end
end

f = @pls_ycgosu_v2;
[Xscores, Yscores, Xweights, Yweights, beta, ~, numFeats] = pls_ycgosu_v2(X, Y, 'numcomps', ncomp);
stats.numFeaturesExtracted = numFeats;
stats.Algorithm_name = 'NIPALS';

if ncomp == 0;ncomp = numFeats;end 

tic;
% if simpls
%     not yet on SIMPLS algorithm...!
%     error('SIMPLS algorithm is NOT POSSIBLE yet.....!!!')
%     fprintf('CACULATING BETA with SIMPLS.........!\n')
%     f = @plsregress;
%     [Xweights, Yweights, Xscores, Yscores, beta, ~,~,~] = f(X,Y,ncomp);
%     stats.descript = 'SIMPLS';
% elseif nipals
%     fprintf('CACULATING BETA with NIPALS.........!\n')
%     f = @pls_ycgosu_v2;
%     [Xscores, Yscores, Xweights, Yweights, beta, ~, numFeats] = f(X, Y, 'numcomps', ncomp);
%     % 
%     stats.numFeaturesExtracted = numFeats;
%     stats.descript = 'NIPALS';
%end
toc2hms = sec2hms(toc);
disp(toc2hms)

stats.intercept = beta(1, :);
stats.beta = beta(2:end, :);
yfit_temp = zeros(size(Y));
mse_temp = zeros(size(Y));


teIdx=cvpart.teIdx; trIdx=cvpart.trIdx;

for fold = 1:max(nfolds)

    fprintf([repmat('-', 1, 20), '%03d folds out of %03d', repmat('-', 1, 20), '\n'], fold, max(nfolds))

    xtrain = X(trIdx{fold}, :);
    ytrain = Y(trIdx{fold}, :);
    xtest = X(teIdx{fold}, :);
    ytest = Y(teIdx{fold}, :);

    testidx = find(teIdx{fold});

    [~, ~, ~, ~, beta, ~, ~] = f(xtrain, ytrain, 'numcomps', ncomp);

    yhat = [ones(size(xtest,1), 1) xtest] * beta;
    err = ytest - yhat;

    yfit_temp(testidx, :) = yhat;
    mse_temp(testidx, :) = err;

end

stats.Y = Y;
stats.yfit = yfit_temp;
stats.mse = sum((mse_temp.^2) / size(mse_temp, 1));
[r, ~] = corr(yfit_temp, Y);
stats.pred_outcome_r = diag(r)';


end



