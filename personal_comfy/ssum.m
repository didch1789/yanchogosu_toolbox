function ssumout = ssum(dat)
numsum = numel(size(dat));
ssumout = sum(x, 1:numsum);
end

function [mminout, mminoutidx] = mmin(dat)
mminout = min(dat, [], 'all');
mminoutidx = (min(dat, [], 'all') == dat);
end

function [mmaxout, mmaxoutidx] = mmax(dat)
mmaxout = max(dat, [], 'all');
mmaxoutidx = (max(dat, [], 'all') == dat);
end