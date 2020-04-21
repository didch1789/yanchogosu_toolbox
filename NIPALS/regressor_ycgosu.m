function regressors = regressor_ycgosu(TR, filter, single_vox)

if size(single_vox, 2) == 1
    fy = fft(single_vox - mean(single_vox), [], 1);
    t = (0:TR:(size(single_vox, 1)-1)*TR);

    fy = fy(1:length(fy)/2);
    freq = 1/TR * (0:length(t)/2-1)/length(t);

    % plot(freq_domain, real(fy))

    idx = zeros(numel(freq), 1);
    for i = 1:numel(freq)
        if freq(i) < filter(1)
            idx(i) = 1;

        elseif freq(i) > filter(2)
            idx(i) = 1;
        end
    end
    
    idx(1) = 0;
    
    for i = 1:numel(idx)
        if idx(i) == 1
            new_idx = zeros(1, numel(freq));
            new_idx(i) = 1;
            regressor(i, :) = real(ifft(fy' .* new_idx, numel(t)));
        end
    end

    regressors = regressor(logical(idx), :)';
else
    error('input should be a single voxel(column wise)')
end
