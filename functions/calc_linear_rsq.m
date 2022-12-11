function [out, normvec, dim1prj] = calc_linear_rsq(A, Amdl, varargin)
% Inputs
%   A (#obs x #feature). going to be projection in 1 dimension
%       according to the 'distance' in varargin.
%   Amdl (#obs x #feature). Regressor for 1 dim prj. do not include intercept.
% varargin
%   {'distance', str}. see Distance in 'pdist'. Default : euclidean
%
% Outputs
%   out.rsq     : R-squared.
%   out.rsq_adj : adjusted R-squared.
%   normvec     : projection vector in 1dim.

dist_metric = 'euclidean';
for i = 1:numel(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'distance'}
                dist_metric = varargin{i+1};
        end
    end
end

D = squareform(pdist(A, dist_metric));
[R1, ~] = find(D == max(pdist(A, dist_metric)));
idx1 = R1(1); idx2 = R1(2);
dim1vec     = A(idx1, :) - A(idx2, :);
normvec     = dim1vec ./ norm(dim1vec);
dim1prj     = A * normvec';
mdl         = fitlm(Amdl, dim1prj);
out.rsq     = mdl.Rsquared.Ordinary;
out.rsq_adj = mdl.Rsquared.Adjusted;



end