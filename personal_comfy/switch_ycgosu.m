function out = switch_ycgosu(A, idx1, idx2)
% Inputs:
%   - A: target Matrix
%   - idx1: idx of target (should be logical idx same size)
%   - idx2: idx of target (should be logical idx)
%       A, idx1 and 2 should have a same size. and sum of idx1 and idx2
%       should be same
%
% Outputs:
%   - out: changed matrix
%
% Example:
%   A = [2     3     4
%        4     5     2
%        4     2     3
%        1     3     4
%        1     2     5];
%   idx1 = [1 0 1
%           0 0 0
%           1 0 0
%           0 0 0
%           0 0 0];
%   idx2 = [0 1 0
%           0 1 1
%           0 0 0
%           0 0 0
%           0 0 0];
%   out = switch_ycgosu(A, logical(idx1), logical(idx2))
%   out = [     3     2     2
%               4     4     4
%               5     2     3
%               1     3     4
%               1     2     5;

if ~isequal(size(A), size(idx1), size(idx2))
    error('Wrong input size!')
end

idx1dim = numel(size(idx1));
idx2dim = numel(size(idx2));

if ((sum(idx1 == 0, 1:idx1dim) + sum(idx1 ~= 0, 1:idx1dim) == numel(idx1)) && ...
        (sum(idx2 == 0, 1:idx1dim) + sum(idx2 ~= 0, 1:idx2dim) == numel(idx2)))
else
    error('Two input indeces should be only contains 1 or 0!')
end

idx1 = logical(idx1);idx2 = logical(idx2);

val1 = A(idx1);
val2 = A(idx2);

A(idx1) = val2;
A(idx2) = val1;

out = A;
end






