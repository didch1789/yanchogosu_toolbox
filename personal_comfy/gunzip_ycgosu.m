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
    

'/media/das/7T/data/PPD_200619_Hahnz/NII/RUN1/gunzip_ycgosu.m'