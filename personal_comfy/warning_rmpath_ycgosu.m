function warning_rmpath_ycgosu(saveornot)
% copy your warning message and just run the function!
% saveornot 

w_temp = clipboard('paste');
w = split(w_temp, newline); 

for i = 1:numel(w)
    filesepidx = regexp(w{i}, filesep);
    if ~isempty(filesepidx)
        rmpath(deblank(w{i}(filesepidx(1):end)))
        disp([deblank(w{i}(filesepidx(1):end)), ' deleted from the path...!'])
    end
end

if saveornot == 1
    savepath
end