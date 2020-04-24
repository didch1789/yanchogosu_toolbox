function warning_rmpath_ycgosu(saveornot)
% copy your warning message and just run the function!
% saveornot 

replace_copy_ycgosu('Warning: Name is nonexistent or not a directory:', sprintf(''))
replace_copy_ycgosu(newline, ' ')
w = clipboard('paste');

w2 = split(w, ' ');

for i = 1:numel(w2)
    if ~isempty(w2{i})
        rmpath(w2{i})
        disp([w2{2}, 'deleted from the path...!'])
    end
end

if saveornot == 1
    savepath
end


 