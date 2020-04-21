function tobeprint = sprintf_ycgosu(words, mark, varargin)

% input
    % words: words to be printed. ('char')
    % mark: mark to be printed with. ('char')
    
    % varargin
        % 'num': length of total string. default: 90
        % 'loc': 'left', 'center', or 'right'. default: 'center'
        % 'off': remove '\n' in the last sentence. default: '\n' in the last sentence
        % last varargin should be words that you want to substitute with %
    
% output
    % tobeprint: string that contains words and mark.

% Example
    % out = sprint_ycgosu('print_ex_%d', '=', 'num', 60, 30)
    % out = '=======================print_ex_30==========================\n'
    
if ~ischar(mark),error('mark should be a character-!'),end

num = 90;
loc = 'center';
n = 1;

if ~isempty(varargin)
    for i = 1:numel(varargin)
        if ischar(varargin{i})
            switch varargin{i}
                case {'num'}
                    num = varargin{i+1};
                case {'loc'}
                    loc = varargin{i+1};
                case {'off'}
                    n = 0;
            end
        end
    end
end

if contains(words, '%')
   numof = count(words, '%');
   words =  sprintf(words, varargin{end-(numof-1):end});
end

str_template = repmat(mark, 1, num);
word_len = length(words);

switch loc
    case {'center'}
        startidx = floor((num - word_len) / 2);
        str_template(startidx:startidx+word_len-1) = words;
    case {'left'}
        str_template(1:word_len) = words;
    case {'right'}
        str_template(num-(word_len-1):end) = words;
end

if n == 1
    tobeprint = [str_template, '\n'];
elseif n == 0
    tobeprint = str_template;
end

end