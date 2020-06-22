function gunzip_ycgosu(deleteornot)
% do gunzip in current directory
% input:
%   should be logical(1 or 0) and it removes .gz file if logical is
%   true.

A = sort_ycgosu(fullfile(pwd,  '*gz'));

for i = 1:numel(A)
    gunzip(A)
    if deleteornot
        delete(A{i});
    end
end
   