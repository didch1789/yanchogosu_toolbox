function out = regress_CV(X, y, CVS, varargin)
% Input:
%     X : #obs X # feats (It doesn't care about Intercept so it should be
%     added at the input stage).
%     y : #obx X # num 
%     CVS : matlab cvpartition object (e.g., CVS = cvpartition(#obs,'KFold', 5) )
% OutPut:
%   out :
%       .beta : regression coefficients using full data.
%       .corr : correlation
%       .rsq  : r-squared.
%       .arsq : adjusted r-squared.

do_verb = false;
for i = 1:numel(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'verbose'
                do_verb = true;
        end
    end
end

if do_verb
    fprintf('Doing %02d Folds\n', CVS.NumTestSets)
end

[N, p] = size(X);
if any(all(X == 1))
    p = p - 1;
end
out.beta = pinv(X) * y;
yfit_all = X*pinv(X)*y;
out.rsq = (1 - sum((y-yfit_all).^2) ./ sum((y-mean(y)).^2));
out.corr = corr(y, yfit_all);


yfit = NaN(size(y));
for cv = 1:CVS.NumTestSets
    trX = X(training(CVS, cv), :);
    trY = y(training(CVS, cv), :);
    teX = X(test(CVS, cv), :);
    trBeta = pinv(trX) * trY;
    yfit(test(CVS, cv), :) = teX * trBeta;
end

out.rsq_cv =  1 - (sum((y-yfit).^2) ./ sum((y-mean(y)).^2));
out.arsq_cv = 1 - ((1 - out.rsq_cv) * (N - 1) / (N - p - 1));
out.corr_cv = corr(y, yfit);



end