function C = dct_ycgosu(TR, filter, k)
% Creates basis functions for Discrete Cosine Transform.
%
% TR - repetition time
% k - data to identify whole time series (rows: voxels, columns: time series)
% filter - [x y] 
%   x = high pass threshold
%   y = low pass threshold
%__________________________________________________________________________ 
 
% initialise
%--------------------------------------------------------------------------

N = size(k, 2);
C = ones(N,1);
dt = 0:TR:(N-1)*TR;
freq = (1/TR) * (1:length(dt)/2)/length(dt);

for i = 1:numel(freq)-1
    if freq(i) < filter(1)
        C(:,end + 1) = sin(2*pi*freq(i)*dt);
        C(:,end + 1) = cos(2*pi*freq(i)*dt);
    elseif freq(i) > filter(2)
        C(:,end + 1) = sin(2*pi*freq(i)*dt);
        C(:,end + 1) = cos(2*pi*freq(i)*dt);
    end
end

C(:, 1) = [];

end