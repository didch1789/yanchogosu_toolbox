function out = region_insert_ycgosu(in, vals, wh_type)

if nargin < 3
    error('input should be region obj, vals[numel(region) X 1], type')
end

numregs = numel(in);
numvals = numel(vals);

if numvals ~= numregs
    error('number of first and second input should be same!')
end

if wh_type == 'Z'
    tt = 1;
elseif wh_type == 'dat'
    tt = 2;
elseif wh_type == 'val'
    tt = 3;
elseif wh_type == 'all'
    tt = 4;
end

types_temp = {'Z', 'dat', 'val', {'Z', 'dat', 'val'}};
types = types_temp{tt};

for i = 1:numregs
    if ~iscell(types)
        in(:, i).(types) = repmat(vals(i), size(in(:, i).(types)));
    else
        for j = 1:numel(types)
            in(:, i).(types{j}) = repmat(vals(i), size(in(:, i).(types{j})));
        end
    end
end

out = in;
end

