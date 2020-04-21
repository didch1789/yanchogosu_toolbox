function [mmaxout, mmaxoutidx] = mmax(dat)
mmaxout = max(dat, [], 'all');
mmaxoutidx = (max(dat, [], 'all') == dat);
end