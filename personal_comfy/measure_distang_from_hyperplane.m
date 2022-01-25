function [distout, residout, angleout] = measure_distang_from_hyperplane(plane, points, varargin)
% plane should be n x m matrix, which spands column space in n dimesion.
%   (plane is assumed to pass the center origin).
% points should be k x n matrix, which has k data points in n dimension.
% Same using GLM residuals.
% varargin:
%   'addintrcpt' : add ones, which has mean centering effect. default is false

addint = false;
for i = 1:numel(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'addint'
                addint = true;
        end
    end         
end

if addint
    plane = [ones(size(plane, 1), 1) plane];
end

points_hat = plane * pinv(plane) * points';
residout = points_hat;
resids = (points' - points_hat);
distout = sqrt(sum(resids.^2));
angleout = NaN(1, size(resids, 2));
for i = 1:size(resids, 2)
    angleout(:, i) = corr(points(i, :)', points_hat(:, i));
end

% Nullspace = null(plane');
% vectors = Nullspace * pinv(Nullspace) * points';
% distout = sqrt(sum(vectors.^2));


end

