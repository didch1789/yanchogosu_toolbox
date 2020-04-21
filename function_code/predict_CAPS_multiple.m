function pred_model = predict_CAPS_multiple(obj, varargin)

% ..
%    ---------------------------------------------------------------------
%    Defaults
%    ---------------------------------------------------------------------
% ..

nfolds = 5;
tsxval_hvblock = 0;
tsxval_rolling = 0;
tsxval_rolling_CAPS = 0;
%error_type = 'mcr';            % mcr, mse: misclassification rate or mean sq. error
algorithm_name = 'cv_pcr';      % name of m-file defining training/test function
useparallel = 'always';         % Use parallel processing, if available, for bootstrapping. Currently no parallel for xval. always = do it, anything else, don't use parallel
bootweights = 0;                % bootstrap voxel weights
verbose = 1;
doMultiClass = 0;               % option to run multiclass SVM
cv_save = 0;


if length(unique(obj.Y)) == 2
    error_type = 'mcr';
else
    error_type = 'mse';
end

if size(obj.dat, 2) ~= size(obj.Y, 1)
    if ~strcmp('MultiClass',varargin) %Added by LC 2/26/13 for new svm multiclass functionality
        error('obj.dat must be [Predictors x observations] and obj.Y must be [observations x 1]');
    end
end

% Mandatory conditioning of data object
% ---------------------------------------------------------------------

%obj = remove_empty(obj);
%11/27/12: Luke Chang: This seems to be
%affecting the ability to use orthviews.  This should be either removed or
%a 'replace_empty(obj) needs to be added somewhere. % TOR: orthviews should
%be able to handle this by doing replace_empty there if needed

% force double to avoid various problems
obj.dat = double(obj.dat);
obj.Y = double(obj.Y);

% ---------------------------------------------------------------------
% Split inputs into those that control crossval and those that
% belong to the algorithm (can be different for different
% algorithms).
% ---------------------------------------------------------------------

predfun_inputs = {};
bootfun_inputs = {}; % for setting number of bootstrap samples

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            
            case {'noparallel'}
                useparallel = 'never';
                
            case {'nfolds', 'error_type', 'algorithm_name', 'useparallel', 'verbose'}
                str = [varargin{i} ' = varargin{i + 1};'];
                eval(str)
                varargin{i} = [];
                varargin{i + 1} = [];
                
            case {'cv_pcr', 'cv_multregress', 'cv_univregress', 'cv_svr', 'cv_lassopcr', 'cv_svm','cv_lassopcrmatlab'}
                algorithm_name = varargin{i};
                
            case {'bootweights'}
                bootweights = 1; varargin{i} = [];
                
            case 'bootsamples'
                bootfun_inputs{end+1} = 'bootsamples';
                bootfun_inputs{end+1} = varargin{i+1};
                bootweights = 1;
                varargin{i} = [];
                varargin{i+1} = [];
                
            case 'savebootweights'
                bootfun_inputs{end+1} = 'savebootweights';
                bootweights = 1;
                varargin{i} = [];
                
            case 'MultiClass'
                doMultiClass = 1;
                predfun_inputs{end + 1} = varargin{i}; %PK - pass as input to cv_svm
                
            case 'rolling'
                tsxval_rolling = 1;
                inval = varargin{i + 1};
                h = inval(1);
                v = inval(2);
                g = inval(3);
                varargin{i} = [];
                varargin{i + 1} = [];
                
            case 'rolling_CAPS'
                tsxval_rolling_CAPS = 1;
                inval = varargin{i + 1};
                h = inval(1);
                v = inval(2);
                g = inval(3);
                cv_set = inval(4);
                cv_method = inval(5);
                varargin{i} = [];
                varargin{i + 1} = [];
                
            case 'hvblock'
                tsxval_hvblock = 1;
                inval = varargin{i + 1};
                h = inval(1);
                v = inval(2);
                varargin{i} = [];
                varargin{i + 1} = [];
                
            case 'savecv'
                cv_save = 1;
                
            otherwise
                % use all other inputs as optional arguments to specific
                % functions, to be interpreted by them
                predfun_inputs{end + 1} = varargin{i};
                
                if (i+1) <= length(varargin) && ~ischar(varargin{i + 1})
                    predfun_inputs{end + 1} = varargin{i + 1};
                end
                
        end
    end
end

% For CPM
cpm_pos = false;
cpm_neg = false;
cpm_thresh = 0.01;
cpm_spearman = false;
cpm_robust = false;
cpm_partialcorr = false;
partialvec = [];
cpm_polynom = 1;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'polarity'
                if strcmp(varargin{i+1}, 'pos')
                    cpm_pos = true;
                elseif strcmp(varargin{i+1}, 'neg')
                    cpm_neg = true;
                elseif strcmp(varargin{i+1}, 'both')
                    cpm_pos = true;
                    cpm_neg = true;
                end
            case 'threshnum'
                cpm_thresh = varargin{i+1};
            case 'spearman'
                cpm_spearman = true;
            case 'robust'
                cpm_robust = true;
            case 'partialcorr'
                cpm_partialcorr = true;
                partialvec = varargin{i+1};
            case 'polynorm'
                cpm_polynorm = varargin{i+1};
        end
    end
end


% Set up parallel processing input to bootstrap function
switch useparallel
    case 'always'
        bootfun_inputs{end+1} = 'parallel';
    otherwise
        % no parallel
end



% This is not used.
% opt = statset('crossval');
% opt.UseParallel = useparallel;

% ---------------------------------------------------------------------
% stratified partition, or custom holdout set
% ---------------------------------------------------------------------

% NOTE: COULD REPLACE ALL THIS WITH cvpart = stratified_holdout_set(Y, varargin)

if numel(nfolds) > 1 % a vector; assume integers for holdout sets
    % Custom holdout set
    fold_indicator = nfolds;
    u = unique(fold_indicator);
    nfolds = length(u);
    
    cvpart = struct('NumTestSets', nfolds);
    [trIdx, teIdx] = deal(cell(1, nfolds));
    
    for i = 1:length(u)
        teIdx{i} = fold_indicator == u(i);
        trIdx{i} = ~teIdx{i};
    end
    
elseif nfolds == 1 % special for 1 fold: use all obs for train and test; good for bootstrapping weights
    [trIdx, teIdx] = deal(cell(1, nfolds));
    trIdx={ones(size(obj.Y))};
    teIdx={ones(size(obj.Y))};
    cvpart = struct('NumTestSets', length(trIdx)); %LC: 4/3/14: added this as it was missing
    
elseif tsxval_hvblock == 1 %special case for timeseries CV using HVBlock : added LC: 11/28/13
    [trIdx, teIdx] = tscv(length(obj.Y), 'hvblock',[h,v]);
    cvpart = struct('NumTestSets', length(trIdx));
    
elseif tsxval_rolling == 1 %special case for timeseries CV using HVBlock : added LC: 12/16/13
    [trIdx, teIdx] = tscv(length(obj.Y), 'rolling',[h,v,g]);
    cvpart = struct('NumTestSets', length(trIdx));
    
    
elseif tsxval_rolling_CAPS == 1 %special case for timeseries CV using HVBlock : added JJ: 07/16/17
    [trIdx, teIdx] = tscv(length(obj.Y)/cv_set, 'rolling',[h,v,g]);
    if cv_method == 1
        for cv_i = 1:numel(trIdx)
            trIdx{cv_i} = repmat(trIdx{cv_i}, 1, cv_set)';
            trIdx{cv_i} = trIdx{cv_i}(:);
        end
        for cv_i = 1:numel(teIdx)
            teIdx{cv_i} = repmat(teIdx{cv_i}, 1, cv_set)';
            teIdx{cv_i} = teIdx{cv_i}(:);
        end
    elseif cv_method == 2
        for cv_i = 1:numel(trIdx)
            trIdx{cv_i} = repmat(trIdx{cv_i}, cv_set, 1);
        end
        for cv_i = 1:numel(teIdx)
            teIdx{cv_i} = repmat(teIdx{cv_i}, cv_set, 1);
        end
    end
    cvpart = struct('NumTestSets', length(trIdx));
    
else
    
    [trIdx, teIdx] = deal(cell(1, nfolds));
    
    % Classification: Stratified holdout set
    if doMultiClass
        
        %cvpart = cvpartition(obj.Y,'k',nfolds); %changed by LC 2/26/13 to account for multiclass svm
        cvpart = cvpartition(length(obj.Y),'k',nfolds);
        
        for i = 1:cvpart.NumTestSets
            
            trIdx{i} = cvpart.training(i);
            teIdx{i} = cvpart.test(i);
            
        end
        
    else
        % Regression: custom stratification
        % cvpartition object will not stratify continuous values
        % do our own
        
        cvpart = struct('NumTestSets', nfolds);
        
        [ys, indx] = sort(obj.Y);
        for k = 1:nfolds
            
            wh_holdout = indx(k:nfolds:end);
            if isempty(wh_holdout), error('Holdout set construction error: Dataset too small?'); end
            
            teIdx{k} = false(size(obj.Y));
            teIdx{k}(wh_holdout) = true;
            
            trIdx{k} = ~teIdx{k};
        end
    end
end

if verbose
    fprintf('Cross-validated prediction with algorithm %s, %3.0f folds\n', algorithm_name, nfolds)
end

% special for 1 fold: use all obs for train and test; good for
% bootstrapping weights
if nfolds == 1
    %trIdx{1} = teIdx{1}; %LC: this isn't actualy using all of the weights.
    % 3/23/13 Wani Woo added the following three lines to get cv_assignment. See Programmers' note for the detail.
    cv_assignment = double(teIdx{1});
else
    cv_assignment = sum((cat(2, teIdx{:}).*repmat((1:size(cat(2, teIdx{:}), 2)), size(cat(2, teIdx{:}), 1), 1))')';
end

% ---------------------------------------------------------------------
% build function handle, given custom inputs to function
% ---------------------------------------------------------------------

% 3/23/13 Wani Woo added cv_assignment in funstr to feed this variable to cv_lassopcr
funstr = ['@(xtrain, ytrain, xtest, cv_assignment) ' algorithm_name '(xtrain, ytrain, xtest, cv_assignment, predfun_inputs{:})'];

eval(['funhan = ' funstr ';'])


% ---------------------------------------------------------------------
% Cross-validated prediction and error
% ---------------------------------------------------------------------

% Use this loop instead of the crossval function
% it's not complicated and allows more control.
% yfit = crossval(funhan, obj.dat',obj.Y, 'Partition', cvpart, 'Options', opt);


t1 = clock;

% Fit on all data - apparent loss, optional outputs


xtrain = obj.dat';
xtest = obj.dat';
ytrain = obj.Y;



switch algorithm_name
    
    case 'cv_lassopcr_numc'
        
        %%%%%%%%%%%%%%%%%ALGORITHM START%%%%%%%%%%%%%%%%%%%%%%%%
        doSkip = 0; %skip predict and output null \\
        dopcr = 1;
        
        wh = find(strcmp(varargin, 'nopcr'));
        if ~isempty(wh), dopcr = 0; end
        
        if dopcr
            [pc,~,~] = svd(scale(xtrain,1)', 'econ'); % replace princomp with SVD on transpose to reduce running time.
            pc(:,end) = [];
            
            % Choose number of components to save [optional]
            wh = find(strcmp(varargin, 'numcomponents'));
            if ~isempty(wh) && length(varargin) >= wh + 1
                
                numc = varargin{wh + 1};
                
                if numc <= size(pc, 2)
                    % it is correct.
                else
                    disp('WARNING!! Number of components requested is more than unique components in training data. It`ll be fixed.');
                    numc(~[numc <= size(pc, 2)]) = size(pc, 2);
                end
                pc = pc(:, 1:numc);
            end
            
            sc = xtrain * pc;
            
        else
            sc = xtrain;
        end
        
        
        numcomps = rank(sc, .001);
        
        
        if dopcr
            out = lasso_rocha_CAPS(ytrain, sc(:, 1:numcomps));
        else
            out = lasso_rocha_CAPS(ytrain, sc);
        end
        
        n_in_model = sum(out.beta ~= 0, 2); % lambda rows * varaible columns
        
        lassonum = sum(out.beta(end,:) ~= 0); % a number, or an array of number of selection
        
        wh = find(strcmp(varargin, 'lasso_num'));
        if ~isempty(wh) && length(varargin) >= wh + 1
            lassonum = varargin{wh + 1};
            if lassonum <= sum(out.beta(end,:) ~= 0)
                % it is correct.
            else
                disp('You asked for more components/vars than there are in the dataset. It`ll be fixed to max num.')
                lassonum(~[lassonum <= sum(out.beta(end,:) ~= 0)]) = sum(out.beta(end,:) ~= 0);
            end
        elseif ~isempty(wh)
            error('Follow lasso_num input with number of vars to keep');
        end
        
        pred_model = cell(1, numel(lassonum));
        
        for l_i = 1:numel(lassonum)
            
            n = find(n_in_model == lassonum(l_i), 1, 'last'); % find the lambda of least shrinkage
            % sometimes this may not exist, depending on sampling res of b (100 steps by default...)
            % if empty, find the next closest one by relaxing lambda a bit.
            indx = 1;
            while isempty(n)
                n = find(n_in_model == lassonum(l_i)+indx, 1, 'last');
                indx = indx + 1;
            end
            
            wh_beta = logical([[out.intercept(n) ~= 0]; [out.beta(n,:) ~= 0]']);
            btemp = pinv([ones(size(ytrain)) sc(:, wh_beta(2:end))]) * ytrain;
            b = zeros(size(wh_beta));
            b(wh_beta) = btemp;
            
            out_cv = [];
            pred_model{l_i}.optout{3} = out_cv;
            
            if dopcr
                vox_weights = pc(:, 1:numcomps) * b(2:end);
            else
                vox_weights = b(2:end);
            end
            
            intercept = b(1);
            
            pred_model{l_i}.optout{1} = vox_weights;
            pred_model{l_i}.optout{2} = intercept;
            %pred_model{l_i}.stats.yfit = intercept + xtest * vox_weights;
        end
        %%%%%%%%%%%%%%%%%ALGORITHM END%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
    case 'cv_pcr_numc'
        
        %%%%%%%%%%%%%%%%%ALGORITHM START%%%%%%%%%%%%%%%%%%%%%%%%
        [pc,~,~] = svd(scale(xtrain,1)', 'econ'); % replace princomp with SVD on transpose to reduce running time.
        pc(:,end) = [];                % remove the last component, which is close to zero.
        
        
        % Choose number of components to save [optional]
        wh = find(strcmp(varargin, 'numcomponents'));
        if ~isempty(wh) && length(varargin) >= wh + 1
            
            numc = varargin{wh + 1};
            
            if numc <= size(pc, 2)
                % it is correct.
            else
                disp('WARNING!! Number of components requested is more than unique components in training data. It`ll be fixed.');
                numc(~[numc <= size(pc, 2)]) = size(pc, 2);
            end
        end
        
        pred_model = cell(1, numel(numc));
        
        for n_i = 1:numel(numc)
            
            pc_refine = pc(:, 1:numc(n_i));
            sc = xtrain * pc_refine;
            
            numcomps = rank(sc);
            
            % 3/8/13: TW:  edited to use numcomps, because sc is not always full rank during bootstrapping
            X = [ones(size(sc, 1), 1) sc(:, 1:numcomps)];
            
            if rank(X) <= size(sc, 1)
                b = pinv(X) * ytrain; % use pinv to stabilize; not full rank
            else
                b = inv(X'*X)*X'*ytrain;
            end
            
            vox_weights = pc_refine(:, 1:numcomps) * b(2:end);
            
            intercept = b(1);
            
            pred_model{n_i}.optout{1} = vox_weights;
            pred_model{n_i}.optout{2} = intercept;
            %pred_model{n_i}.stats.yfit = intercept + xtest * vox_weights;
            
        end
        %%%%%%%%%%%%%%%%%ALGORITHM END%%%%%%%%%%%%%%%%%%%%%%%%
        
    case 'cv_cpm_numc'
        
        %%%%%%%%%%%%%%%%%ALGORITHM START%%%%%%%%%%%%%%%%%%%%%%%%
        no_observation = size(xtrain, 1);
        no_node = size(xtrain, 2);
        if cpm_partialcorr && size(partialvec, 1) ~= no_node
            disp('controlling variable dimension dismatch with original data. rows should be observations.');
        end
        if cpm_pos && cpm_neg && cpm_polynorm ~= 1
            disp('with both pos and neg fitting, polyfit 2 or more order is banned. We`ll use 1 order multiple linear regression.');
        end
        
        % correlate all edges with behavior
        
        
        if ~cpm_robust
            if cpm_partialcorr && ~cpm_spearman
                [r_mat, p_mat] = partialcorr(xtrain, ytrain, partialvec);
            elseif cpm_partialcorr && cpm_spearman
                [r_mat, p_mat] = partialcorr(xtrain, ytrain, partialvec, 'type', 'spearman');
            elseif ~cpm_partialcorr && ~cpm_spearman
                [r_mat, p_mat] = corr(xtrain, ytrain);
            elseif ~cpm_partialcorr && cpm_spearman
                [r_mat, p_mat] = corr(xtrain, ytrain, 'type', 'spearman');
            end
        elseif cpm_robust
            r_mat = zeros(1, no_node);
            p_mat = zeros(1, no_node);
            for node_i = 1:no_node
                [~, cpm_stats] = robustfit(xtrain(:,node_i), ytrain);
                cur_t = cpm_stats.t(2);
                r_mat(node_i) = sign(cur_t) * sqrt(cur_t^2/(no_observation-1-2+cur_t^2));
                p_mat(node_i) = 2*tcdf(cur_t, no_observation-1-2);
            end
        end
        
        % set threshold and define masks
        th_range = 1:numel(cpm_thresh);
        for t_i = th_range
            issigedge = true;
            if cpm_pos
                pos_edges = find(r_mat > 0 & p_mat < cpm_thresh(t_i));
                if numel(pos_edges) == 0; issigedge = false; end
                wh_pos = r_mat > 0 & p_mat < cpm_thresh(t_i);
            end
            if cpm_neg
                neg_edges = find(r_mat < 0 & p_mat < cpm_thresh(t_i));
                if numel(neg_edges) == 0; issigedge = false; end
                wh_neg = r_mat < 0 & p_mat < cpm_thresh(t_i);
            end
            
            if issigedge
                
                if cpm_pos; fprintf('Threshold %.10f : %.8d/%.8d positive weight edges were detected.\n', cpm_thresh(t_i), numel(pos_edges), no_node); end
                if cpm_neg; fprintf('Threshold %.10f : %.8d/%.8d negative weight edges were detected.\n', cpm_thresh(t_i), numel(neg_edges), no_node); end
                
                % get sum of all edges in TRAIN subs
                if cpm_pos; xtrain_sumpos = sum(xtrain(:, pos_edges), 2); end
                if cpm_neg; xtrain_sumneg = sum(xtrain(:, neg_edges), 2); end
                
                % build model on TRAIN subs
                %     if cpm_pos && cpm_neg % both pos and neg weight
                %         fit_coef = regress(ytrain, [xtrain_sumpos, xtrain_sumneg, ones(no_observation, 1)]); % multiple linear regression
                %         xtest_sumpos = sum(xtest(:, pos_edges), 2);
                %         xtest_sumneg = sum(xtest(:, neg_edges), 2);
                %         yfit = fit_coef(1)*xtest_sumpos + fit_coef(2)*xtest_sumneg + fit_coef(3); % get prediction
                %     elseif cpm_pos && ~cpm_neg % pos weight
                %         fit_coef = polyfit(xtrain_sumpos, ytrain, cpm_polynorm);
                %         xtest_sumpos = sum(xtest(:, pos_edges), 2); % run model on TEST sub
                %         yfit = polyval(fit_coef, xtest_sumpos); % get prediction
                %     elseif ~cpm_pos && cpm_neg % neg weight
                %         fit_coef = polyfit(xtrain_sumneg, ytrain, cpm_polynorm);
                %         xtest_sumneg = sum(xtest(:, neg_edges), 2); % run model on TEST sub
                %         yfit = polyval(fit_coef, xtest_sumneg); % get prediction
                %     end
                
                if cpm_pos && cpm_neg % both pos and neg weight
                    fit_coef = regress(ytrain, [xtrain_sumpos, xtrain_sumneg, ones(no_observation, 1)]); % multiple linear regression
                    vox_weights = fit_coef(1)*wh_pos + fit_coef(2)*wh_neg;
                    intercept = fit_coef(3);
                elseif cpm_pos && ~cpm_neg % pos weight
                    fit_coef = polyfit(xtrain_sumpos, ytrain, cpm_polynorm);
                    vox_weights = polyval(fit_coef, wh_pos) - fit_coef(end);
                    intercept = fit_coef(end);
                elseif ~cpm_pos && cpm_neg % neg weight
                    fit_coef = polyfit(xtrain_sumneg, ytrain, cpm_polynorm);
                    vox_weights = polyval(fit_coef, wh_neg) - fit_coef(end);
                    intercept = fit_coef(end);
                end
                
                pred_model{t_i}.optout{1} = vox_weights;
                pred_model{t_i}.optout{2} = intercept;
                %pred_model{t_i}.stats.yfit = intercept + xtest * vox_weights;
                
            else
                
                fprintf('no significant weight edges were detected. it`ll cause error.\n')
                th_range(th_range == t_i) = [];
                
            end
        end
        %%%%%%%%%%%%%%%%%ALGORITHM END%%%%%%%%%%%%%%%%%%%%%%%%
        
        
end



if verbose
    %fprintf(1,'\n_______________________________\n')
    [hour, minute, second] = sec2hms(etime(clock,t1));
    fprintf(1,'\nCompleted fit for all data in: %3.0f hours %3.0f min %2.0f secs \n',hour,minute,second);
end

% ---------------------------------------------------------------------
% Cross-validated prediction and error
% ---------------------------------------------------------------------


if nfolds ~= 1

    for i = 1:cvpart.NumTestSets

        % 3/23/13 Wani Woo added the following three lines to get a proper cv_assignment variable for cross-validation.

        teIdx_cv = teIdx(1:end ~= i);
        cv_assignment = sum((cat(2, teIdx_cv{:}).*repmat((1:size(cat(2, teIdx_cv{:}), 2)), size(cat(2, teIdx_cv{:}), 1), 1))')';
        cv_assignment(cv_assignment == 0) = [];


        t2 = clock;


        xtrain = obj.dat(:, trIdx{i})';
        xtest = obj.dat(:, teIdx{i})';
        ytrain = obj.Y(trIdx{i},:);

        switch algorithm_name

            case 'cv_lassopcr_numc'
                %%%%%%%%%%%%%%%%%ALGORITHM START%%%%%%%%%%%%%%%%%%%%%%%%

                doSkip = 0; %skip predict and output null \\
                dopcr = 1;

                wh = find(strcmp(varargin, 'nopcr'));
                if ~isempty(wh), dopcr = 0; end

                if dopcr
                    [pc,~,~] = svd(scale(xtrain,1)', 'econ'); % replace princomp with SVD on transpose to reduce running time.
                    pc(:,end) = [];

                    % Choose number of components to save [optional]
                    wh = find(strcmp(varargin, 'numcomponents'));
                    if ~isempty(wh) && length(varargin) >= wh + 1

                        numc = varargin{wh + 1};

                        if numc > size(pc, 2)
                            disp('WARNING!! Number of components requested is more than unique components in training data.');
                            numc = size(pc, 2);
                        end
                        pc = pc(:, 1:numc);
                    end

                    sc = xtrain * pc;

                else
                    sc = xtrain;
                end


                numcomps = rank(sc, .001);


                if dopcr
                    out = lasso_rocha_CAPS(ytrain, sc(:, 1:numcomps));
                else
                    out = lasso_rocha_CAPS(ytrain, sc);
                end

                n_in_model = sum(out.beta ~= 0, 2); % lambda rows * varaible columns

                lassonum = sum(out.beta(end,:) ~= 0); % a number, or an array of number of selection

                wh = find(strcmp(varargin, 'lasso_num'));
                if ~isempty(wh) && length(varargin) >= wh + 1
                    lassonum = varargin{wh + 1};
                    if lassonum <= sum(out.beta(end,:) ~= 0)
                        % it is correct.
                    else
                        disp('You asked for more components/vars than there are in the dataset. It`ll be fixed to max num.')
                        lassonum(~[lassonum <= sum(out.beta(end,:) ~= 0)]) = sum(out.beta(end,:) ~= 0);
                    end
                elseif ~isempty(wh)
                    error('Follow lasso_num input with number of vars to keep');
                end

                for l_i = 1:numel(lassonum)

                    n = find(n_in_model == lassonum(l_i), 1, 'last'); % find the lambda of least shrinkage
                    % sometimes this may not exist, depending on sampling res of b (100 steps by default...)
                    % if empty, find the next closest one by relaxing lambda a bit.
                    indx = 1;
                    while isempty(n)
                        n = find(n_in_model == lassonum(l_i)+indx, 1, 'last');
                        indx = indx + 1;
                    end

                    wh_beta = logical([[out.intercept(n) ~= 0]; [out.beta(n,:) ~= 0]']);
                    btemp = pinv([ones(size(ytrain)) sc(:, wh_beta(2:end))]) * ytrain;
                    b = zeros(size(wh_beta));
                    b(wh_beta) = btemp;

                    if dopcr
                        vox_weights = pc(:, 1:numcomps) * b(2:end); % vox_weights
                    else
                        vox_weights = b(2:end); % vox_weights
                    end

                    intercept = b(1); % intercept

                    pred_model{l_i}.stats.yfit(teIdx{i},:) = intercept + xtest * vox_weights;

                    if cv_save
                        cv_optout{l_i}{i, :} = vox_weights;
                    end

                end
                %%%%%%%%%%%%%%%%%ALGORITHM END%%%%%%%%%%%%%%%%%%%%%%%%


            case 'cv_pcr_numc'

                %%%%%%%%%%%%%%%%%ALGORITHM START%%%%%%%%%%%%%%%%%%%%%%%%
                [pc,~,~] = svd(scale(xtrain,1)', 'econ'); % replace princomp with SVD on transpose to reduce running time.
                pc(:,end) = [];                % remove the last component, which is close to zero.


                % Choose number of components to save [optional]
                wh = find(strcmp(varargin, 'numcomponents'));
                if ~isempty(wh) && length(varargin) >= wh + 1

                    numc = varargin{wh + 1};

                    if numc <= size(pc, 2)
                        % it is correct.
                    else
                        disp('WARNING!! Number of components requested is more than unique components in training data. It`ll be fixed.');
                        numc(~[numc <= size(pc, 2)]) = size(pc, 2);

                    end
                end

                for n_i = 1:numel(numc)

                    pc_refine = pc(:, 1:numc(n_i));
                    sc = xtrain * pc_refine;

                    numcomps = rank(sc);

                    % 3/8/13: TW:  edited to use numcomps, because sc is not always full rank during bootstrapping
                    X = [ones(size(sc, 1), 1) sc(:, 1:numcomps)];

                    if rank(X) <= size(sc, 1)
                        b = pinv(X) * ytrain; % use pinv to stabilize; not full rank
                    else
                        b = inv(X'*X)*X'*ytrain;
                    end

                    vox_weights = pc_refine(:, 1:numcomps) * b(2:end);

                    intercept = b(1);

                    pred_model{n_i}.stats.yfit(teIdx{i},:) = intercept + xtest * vox_weights;

                    if cv_save
                        cv_optout{n_i}{i, :} = vox_weights;
                    end

                end
                %%%%%%%%%%%%%%%%%ALGORITHM END%%%%%%%%%%%%%%%%%%%%%%%%

            case 'cv_cpm_numc'

                %%%%%%%%%%%%%%%%%ALGORITHM START%%%%%%%%%%%%%%%%%%%%%%%%
                no_observation = size(xtrain, 1);
                no_node = size(xtrain, 2);
                if cpm_partialcorr && size(partialvec, 1) ~= no_node
                    disp('controlling variable dimension dismatch with original data. rows should be observations.');
                end
                if cpm_pos && cpm_neg && cpm_polynorm ~= 1
                    disp('with both pos and neg fitting, polyfit 2 or more order is banned. We`ll use 1 order multiple linear regression.');
                end

                % correlate all edges with behavior


                if ~cpm_robust
                    if cpm_partialcorr && ~cpm_spearman
                        [r_mat, p_mat] = partialcorr(xtrain, ytrain, partialvec);
                    elseif cpm_partialcorr && cpm_spearman
                        [r_mat, p_mat] = partialcorr(xtrain, ytrain, partialvec, 'type', 'spearman');
                    elseif ~cpm_partialcorr && ~cpm_spearman
                        [r_mat, p_mat] = corr(xtrain, ytrain);
                    elseif ~cpm_partialcorr && cpm_spearman
                        [r_mat, p_mat] = corr(xtrain, ytrain, 'type', 'spearman');
                    end
                elseif cpm_robust
                    r_mat = zeros(1, no_node);
                    p_mat = zeros(1, no_node);
                    for node_i = 1:no_node
                        [~, cpm_stats] = robustfit(xtrain(:,node_i), ytrain);
                        cur_t = cpm_stats.t(2);
                        r_mat(node_i) = sign(cur_t) * sqrt(cur_t^2/(no_observation-1-2+cur_t^2));
                        p_mat(node_i) = 2*tcdf(cur_t, no_observation-1-2);
                    end
                end

                % set threshold and define masks
                for t_i = th_range
                    issigedge = true;
                    if cpm_pos
                        pos_edges = find(r_mat > 0 & p_mat < cpm_thresh(t_i));
                        if numel(pos_edges) == 0; issigedge = false; end
                        wh_pos = r_mat > 0 & p_mat < cpm_thresh(t_i);
                    end
                    if cpm_neg
                        neg_edges = find(r_mat < 0 & p_mat < cpm_thresh(t_i));
                        if numel(neg_edges) == 0; issigedge = false; end
                        wh_neg = r_mat < 0 & p_mat < cpm_thresh(t_i);
                    end

                    if issigedge

                        if cpm_pos; fprintf('Threshold %.10f : %.8d/%.8d positive weight edges were detected.\n', cpm_thresh(t_i), numel(pos_edges), no_node); end
                        if cpm_neg; fprintf('Threshold %.10f : %.8d/%.8d negative weight edges were detected.\n', cpm_thresh(t_i), numel(neg_edges), no_node); end

                        % get sum of all edges in TRAIN subs
                        if cpm_pos; xtrain_sumpos = sum(xtrain(:, pos_edges), 2); end
                        if cpm_neg; xtrain_sumneg = sum(xtrain(:, neg_edges), 2); end

                        % build model on TRAIN subs
                        %     if cpm_pos && cpm_neg % both pos and neg weight
                        %         fit_coef = regress(ytrain, [xtrain_sumpos, xtrain_sumneg, ones(no_observation, 1)]); % multiple linear regression
                        %         xtest_sumpos = sum(xtest(:, pos_edges), 2);
                        %         xtest_sumneg = sum(xtest(:, neg_edges), 2);
                        %         yfit = fit_coef(1)*xtest_sumpos + fit_coef(2)*xtest_sumneg + fit_coef(3); % get prediction
                        %     elseif cpm_pos && ~cpm_neg % pos weight
                        %         fit_coef = polyfit(xtrain_sumpos, ytrain, cpm_polynorm);
                        %         xtest_sumpos = sum(xtest(:, pos_edges), 2); % run model on TEST sub
                        %         yfit = polyval(fit_coef, xtest_sumpos); % get prediction
                        %     elseif ~cpm_pos && cpm_neg % neg weight
                        %         fit_coef = polyfit(xtrain_sumneg, ytrain, cpm_polynorm);
                        %         xtest_sumneg = sum(xtest(:, neg_edges), 2); % run model on TEST sub
                        %         yfit = polyval(fit_coef, xtest_sumneg); % get prediction
                        %     end

                        if cpm_pos && cpm_neg % both pos and neg weight
                            fit_coef = regress(ytrain, [xtrain_sumpos, xtrain_sumneg, ones(no_observation, 1)]); % multiple linear regression
                            vox_weights = fit_coef(1)*wh_pos + fit_coef(2)*wh_neg;
                            intercept = fit_coef(3);
                        elseif cpm_pos && ~cpm_neg % pos weight
                            fit_coef = polyfit(xtrain_sumpos, ytrain, cpm_polynorm);
                            vox_weights = polyval(fit_coef, wh_pos) - fit_coef(end);
                            intercept = fit_coef(end);
                        elseif ~cpm_pos && cpm_neg % neg weight
                            fit_coef = polyfit(xtrain_sumneg, ytrain, cpm_polynorm);
                            vox_weights = polyval(fit_coef, wh_neg) - fit_coef(end);
                            intercept = fit_coef(end);
                        end

                        pred_model{t_i}.stats.yfit(teIdx{i},:) = intercept + xtest * vox_weights;

                        if cv_save
                            cv_optout{t_i}{i, :} = vox_weights;
                        end

                    else

                        fprintf('no significant weight edges were detected. it`ll get an error.')
                        th_range(th_range == t_i) = [];

                    end
                end
                %%%%%%%%%%%%%%%%%ALGORITHM END%%%%%%%%%%%%%%%%%%%%%%%%

        end

        % 3/23/13 Wani Woo added cv_assignment.
        if verbose
            [hour, minute, second] = sec2hms(etime(clock,t2));
            fprintf(1,['Fold ' num2str(i) '/' num2str(cvpart.NumTestSets) ' done in: %3.0f hours %3.0f min %2.0f sec\n'],hour,minute,second);
        end


    end
    
end


if verbose
    [hour, minute, second] = sec2hms(etime(clock,t1));
    fprintf(1,'\nTotal Elapsed Time = %3.0f hours %3.0f min %2.0f sec\n',hour, minute, second);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------------------------------------------------------------------
% Get error
% ---------------------------------------------------------------------


if nfolds ~= 1

    switch algorithm_name

        case 'cv_lassopcr_numc'

            for l_i = 1:numel(lassonum)

                switch error_type
                    case {'mcr', 'class_loss'}
                        %err = obj.Y ~= yfit;
                        err = obj.Y ~= round(pred_model{l_i}.stats.yfit); %10/7/12: Luke Chang: this will allow mcr to also be calculated using predicted probabilities

                        %cverr = sum(err) ./ length(err);
                        cverr = sum(err) ./ length(err);

                        phi = corr(obj.Y, pred_model{l_i}.stats.yfit); %10/7/12: Luke Chang: this will calculate phi correlation coefficient between two binary variables

                    case {'mse' 'rmse', 'meanabserr'}
                        err = obj.Y - pred_model{l_i}.stats.yfit;

                        %mse = mean(err' * err); %10/8/12: Luke Chang: I think this is only capturing sum of squared error
                        mse = (err' * err)/length(err);  %This should be correct calculation of mean squared error
                        rmse = sqrt(mse);
                        meanabserr = nanmean(abs(err));  %if you are getting strange error here make sure yo are using matlab default nanmean (e.g., which nanmean)
                        r = corrcoef(obj.Y, pred_model{l_i}.stats.yfit, 'rows', 'pairwise');
                        r = r(1, 2);

                        eval(['cverr = ' error_type ';']);

                    otherwise
                        error('Illegal loss function type')
                end


                % ---------------------------------------------------------------------
                % collect output
                % ---------------------------------------------------------------------

                pred_model{l_i}.cverr = cverr;
                pred_model{l_i}.stats = struct('Y', obj.Y, 'algorithm_name', algorithm_name, ...
                    'function_call', funstr, 'function_handle', funhan, ...
                    'yfit', pred_model{l_i}.stats.yfit, 'err', err, 'error_type', error_type, 'cverr', cverr, ...
                    'nfolds', 'nfolds', 'cvpartition', cvpart);

                pred_model{l_i}.stats.teIdx = teIdx;
                pred_model{l_i}.stats.trIdx = trIdx;

                pred_model{l_i}.stats.other_output = []; %optout;
                pred_model{l_i}.stats.other_output_descrip = 'Other output from algorithm - trained on all data (these depend on algorithm)';
                if cv_save
                    pred_model{l_i}.stats.other_output_cv = cv_optout{l_i}; %cv_optout;
                end
                pred_model{l_i}.stats.other_output_cv_descrip = 'Other output from algorithm - for each CV fold';

                %LC: 4/3/14: Removed this as it is causing errors and seems redundant with line 732.
                % %For SVM reorder distance from hyperplane: 12/6/12: Luke Chang: not very elegant solution
                % if strcmp(algorithm_name,'cv_svm')
                %     stats.other_output_cv{3}=ones(length(stats.Y),1);
                %     for i = 1:cvpart.NumTestSets
                %         stats.other_output{3}(find(teIdx{i}==1)) = stats.other_output_cv{i,2};
                %     end
                % end

                switch error_type
                    case {'mcr', 'class_loss'}

                        pred_model{l_i}.stats.phi = phi;

                        if strcmp(algorithm_name,'cv_lassopcrmatlab')
                            pred_model{l_i}.stats.me = mean(abs(obj.Y-pred_model{l_i}.stats.yfit)); %calculate average distance from probability and outcome - should be more sensitive measure of accuracy than mcr for probabilities
                        end

                    case {'mse' 'rmse', 'meanabserr'}

                        pred_model{l_i}.stats.mse = mse;
                        pred_model{l_i}.stats.rmse = rmse;
                        pred_model{l_i}.stats.meanabserr = meanabserr;
                        pred_model{l_i}.stats.pred_outcome_r = r;
                end


                pred_model{l_i}.stats.weight_obj = [];

                pred_model{l_i} = orderfields(pred_model{l_i}, {'cverr', 'stats', 'optout'});

            end


        case 'cv_pcr_numc'

            for n_i = 1:numel(numc)

                switch error_type
                    case {'mcr', 'class_loss'}
                        %err = obj.Y ~= yfit;
                        err = obj.Y ~= round(pred_model{n_i}.stats.yfit); %10/7/12: Luke Chang: this will allow mcr to also be calculated using predicted probabilities

                        %cverr = sum(err) ./ length(err);
                        cverr = sum(err) ./ length(err);

                        phi = corr(obj.Y, pred_model{n_i}.stats.yfit); %10/7/12: Luke Chang: this will calculate phi correlation coefficient between two binary variables

                    case {'mse' 'rmse', 'meanabserr'}
                        err = obj.Y - pred_model{n_i}.stats.yfit;

                        %mse = mean(err' * err); %10/8/12: Luke Chang: I think this is only capturing sum of squared error
                        mse = (err' * err)/length(err);  %This should be correct calculation of mean squared error
                        rmse = sqrt(mse);
                        meanabserr = nanmean(abs(err));  %if you are getting strange error here make sure yo are using matlab default nanmean (e.g., which nanmean)
                        r = corrcoef(obj.Y, pred_model{n_i}.stats.yfit, 'rows', 'pairwise');
                        r = r(1, 2);

                        eval(['cverr = ' error_type ';']);

                    otherwise
                        error('Illegal loss function type')
                end


                % ---------------------------------------------------------------------
                % collect output
                % ---------------------------------------------------------------------

                pred_model{n_i}.cverr = cverr;
                pred_model{n_i}.stats = struct('Y', obj.Y, 'algorithm_name', algorithm_name, ...
                    'function_call', funstr, 'function_handle', funhan, ...
                    'yfit', pred_model{n_i}.stats.yfit, 'err', err, 'error_type', error_type, 'cverr', cverr, ...
                    'nfolds', 'nfolds', 'cvpartition', cvpart);

                pred_model{n_i}.stats.teIdx = teIdx;
                pred_model{n_i}.stats.trIdx = trIdx;

                pred_model{n_i}.stats.other_output = []; %optout;
                pred_model{n_i}.stats.other_output_descrip = 'Other output from algorithm - trained on all data (these depend on algorithm)';
                if cv_save
                    pred_model{n_i}.stats.other_output_cv = cv_optout{n_i};
                end
                pred_model{n_i}.stats.other_output_cv_descrip = 'Other output from algorithm - for each CV fold';

                %LC: 4/3/14: Removed this as it is causing errors and seems redundant with line 732.
                % %For SVM reorder distance from hyperplane: 12/6/12: Luke Chang: not very elegant solution
                % if strcmp(algorithm_name,'cv_svm')
                %     stats.other_output_cv{3}=ones(length(stats.Y),1);
                %     for i = 1:cvpart.NumTestSets
                %         stats.other_output{3}(find(teIdx{i}==1)) = stats.other_output_cv{i,2};
                %     end
                % end

                switch error_type
                    case {'mcr', 'class_loss'}

                        pred_model{n_i}.stats.phi = phi;

                        if strcmp(algorithm_name,'cv_lassopcrmatlab')
                            pred_model{n_i}.stats.me = mean(abs(obj.Y-pred_model{n_i}.stats.yfit)); %calculate average distance from probability and outcome - should be more sensitive measure of accuracy than mcr for probabilities
                        end

                    case {'mse' 'rmse', 'meanabserr'}

                        pred_model{n_i}.stats.mse = mse;
                        pred_model{n_i}.stats.rmse = rmse;
                        pred_model{n_i}.stats.meanabserr = meanabserr;
                        pred_model{n_i}.stats.pred_outcome_r = r;
                end



                pred_model{n_i}.stats.weight_obj = [];

                pred_model{n_i} = orderfields(pred_model{n_i}, {'cverr', 'stats', 'optout'});
            end

        case 'cv_cpm_numc'

            pred_model = pred_model(th_range);

            for t_i = th_range

                switch error_type
                    case {'mcr', 'class_loss'}
                        %err = obj.Y ~= yfit;
                        err = obj.Y ~= round(pred_model{t_i}.stats.yfit); %10/7/12: Luke Chang: this will allow mcr to also be calculated using predicted probabilities

                        %cverr = sum(err) ./ length(err);
                        cverr = sum(err) ./ length(err);

                        phi = corr(obj.Y, pred_model{t_i}.stats.yfit); %10/7/12: Luke Chang: this will calculate phi correlation coefficient between two binary variables

                    case {'mse' 'rmse', 'meanabserr'}
                        err = obj.Y - pred_model{t_i}.stats.yfit;

                        %mse = mean(err' * err); %10/8/12: Luke Chang: I think this is only capturing sum of squared error
                        mse = (err' * err)/length(err);  %This should be correct calculation of mean squared error
                        rmse = sqrt(mse);
                        meanabserr = nanmean(abs(err));  %if you are getting strange error here make sure yo are using matlab default nanmean (e.g., which nanmean)
                        r = corrcoef(obj.Y, pred_model{t_i}.stats.yfit, 'rows', 'pairwise');
                        r = r(1, 2);

                        eval(['cverr = ' error_type ';']);

                    otherwise
                        error('Illegal loss function type')
                end


                % ---------------------------------------------------------------------
                % collect output
                % ---------------------------------------------------------------------

                pred_model{t_i}.cverr = cverr;
                pred_model{t_i}.stats = struct('Y', obj.Y, 'algorithm_name', algorithm_name, ...
                    'function_call', funstr, 'function_handle', funhan, ...
                    'yfit', pred_model{t_i}.stats.yfit, 'err', err, 'error_type', error_type, 'cverr', cverr, ...
                    'nfolds', 'nfolds', 'cvpartition', cvpart);

                pred_model{t_i}.stats.teIdx = teIdx;
                pred_model{t_i}.stats.trIdx = trIdx;

                pred_model{t_i}.stats.other_output = []; %optout;
                pred_model{t_i}.stats.other_output_descrip = 'Other output from algorithm - trained on all data (these depend on algorithm)';
                if cv_save
                    pred_model{t_i}.stats.other_output_cv = cv_optout{t_i};
                end
                pred_model{t_i}.stats.other_output_cv_descrip = 'Other output from algorithm - for each CV fold';

                %LC: 4/3/14: Removed this as it is causing errors and seems redundant with line 732.
                % %For SVM reorder distance from hyperplane: 12/6/12: Luke Chang: not very elegant solution
                % if strcmp(algorithm_name,'cv_svm')
                %     stats.other_output_cv{3}=ones(length(stats.Y),1);
                %     for i = 1:cvpart.NumTestSets
                %         stats.other_output{3}(find(teIdx{i}==1)) = stats.other_output_cv{i,2};
                %     end
                % end

                switch error_type
                    case {'mcr', 'class_loss'}

                        pred_model{t_i}.stats.phi = phi;

                        if strcmp(algorithm_name,'cv_lassopcrmatlab')
                            pred_model{t_i}.stats.me = mean(abs(obj.Y-pred_model{t_i}.stats.yfit)); %calculate average distance from probability and outcome - should be more sensitive measure of accuracy than mcr for probabilities
                        end

                    case {'mse' 'rmse', 'meanabserr'}

                        pred_model{t_i}.stats.mse = mse;
                        pred_model{t_i}.stats.rmse = rmse;
                        pred_model{t_i}.stats.meanabserr = meanabserr;
                        pred_model{t_i}.stats.pred_outcome_r = r;
                end


                
                pred_model{t_i}.stats.weight_obj = [];
                
                pred_model{t_i} = orderfields(pred_model{t_i}, {'cverr', 'stats', 'optout'});
            end
    end
    
end

end % main function


% ------------ bootstrapping and other functions -------------------------

function [hour, minute, second] = sec2hms(sec)
%SEC2HMS  Convert seconds to hours, minutes and seconds.
%
%   [HOUR, MINUTE, SECOND] = SEC2HMS(SEC) converts the number of seconds in
%   SEC into hours, minutes and seconds.

hour   = fix(sec/3600);      % get number of hours
sec    = sec - 3600*hour;    % remove the hours
minute = fix(sec/60);        % get number of minutes
sec    = sec - 60*minute;    % remove the minutes
second = sec;
end

% ------------ bootstrapping and other functions -------------------------

function rho = nancorr(a,b)
%calculate pearson correlation between a & b for nonnan values
keep = logical(~isnan(a) .* ~isnan(b));
rho = (mean((a(keep)-mean(a(keep))) .* (b(keep)-mean(b(keep)))))/(std(a(keep))*std(b(keep)));
end
