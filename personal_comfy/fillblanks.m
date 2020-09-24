function out = fillblanks(varargin)
% exact same logic with function 'fullfile'
% Difference is that fullfile fills with '/' while
% fillblanks fills with blanks(1)

out = strrep(fullfile(varargin), filesep, blanks(1));

end