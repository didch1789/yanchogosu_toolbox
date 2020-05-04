function [out1, out2] = fileparts_ycgosu(filedir, howmany)
% Input options:
%   "filedir" can be a char, string or cell array.
%   "howmany" defines number of back slash(integer)
%
% outputs:
% "out1" gives a filedir before "howmany"
% "out2" gives a filedir after "howmany" (include stuff in howmany and filename)
%
% fileparts_ycgosu(filedir, 1) is equivalent to fileparts(filedir) expect
% that fileparts seperates the filename and ext.

switch class(filedir)
    case 'cell'
        target = filedir;
        docell = 1;
        dochar = 0;
    case {'string', 'char'}
        if size(filedir, 1) > 1
            disp('Rather put it in cell array!!')
        else
            target = filedir;
            docell = 0;
            dochar = 1;
        end
    otherwise
        disp('Unknown input!')
end

template = split(target, filesep);

if dochar
    out1 = strjoin(template(1:(end-howmany)), filesep);
    out2 = strjoin(template(end-(howmany-1):end), filesep);
elseif docell
    out1_temp = template(:, 1:(end-howmany));
    out2_temp = template(:, end-(howmany-1):end);
    for i = 1:size(out1_temp, 1)
        out1{i, :} = strjoin(out1_temp(i, :), filesep);
        out2{i, :} = strjoin(out2_temp(i, :), filesep);
    end
end
    

end