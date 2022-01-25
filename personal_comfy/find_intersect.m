function [outsect, tol] = find_intersect(plane1, plane2, varargin)
% plane1 and plane2 should be one that shares same dimension and spanning column space of each.
% plane1, plane2 : n x k (n dimension, with k columns)

numC = 3;
for i = 1:numel(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'numC'
                numC = varargin{i+1};
        end
    end
end

prj1 = plane1 * pinv(plane1);
prj2 = plane2 * pinv(plane2);
PRJ = (prj1 + prj2);
A = (PRJ+PRJ')/2;

[eigVec, eigVal] = eig(A);
reigVal = diag(eigVal);
[maxvals, maxidx] = maxk(reigVal, numC);
outsect = eigVec(:, maxidx);
tol     = 2 - maxvals;

end
