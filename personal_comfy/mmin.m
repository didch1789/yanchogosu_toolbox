function [mminout, mminoutidx] = mmin(dat)
mminout = min(dat, [], 'all');
mminoutidx = (min(dat, [], 'all') == dat);
end