classdef spk_data
    properties
        spikes;
        pseudopops;
        dims = 'trial X neuron X time'
        history = {'Processed History:'};
    end
    
    methods
        function obj = binning(obj, binsize)
            rasterdat = obj.spikes;
            if numel(size(rasterdat)) > 2 
                n_trial  = size(rasterdat, 1);
                n_neuron = size(rasterdat, 2);
                n_time   = size(rasterdat, 3);
                Q        = fix(n_time / binsize);
                R        = mod(n_time, binsize);
                binned   = NaN(n_trial, n_neuron, Q);

                for tt = 1:n_trial
                    spktrials   = squeeze(rasterdat(tt, :, :));
                    l = 1;
                    for bb = 1:Q
                        if R ~= 0
                            k = binsize + 1;
                            R = R - 1;
                        else
                            k = binsize;
                        end
                        binned(tt, :, bb) = sum(spktrials(:, l:l+k-1), 2);
                        l = l + k;
                    end
                end
            else
                error('Rather use bin_raster_ycgosu.m');
            end
            obj.spikes = binned;
            obj.history{end+1} = sprintf('binned with size %02d', binsize);
        end

        function obj = smoothing(obj, windowsize)
            obj.spikes = smoothdata(obj.spikes, 3, 'gaussian', windowsize);
            obj.history{end+1} = sprintf('smoothed with size %d', windowsize);
        end

        function obj = condition_avg(obj, condIds, varargin)
            n_trial = size(obj.spikes, 1);
            if n_trial ~= numel(condIds)
                error('Check number of trials')
            end
            
            n_cond = unique(numel(condIds));
            for i = 1:n_cond
                obj.pseudopops{i} = squeeze(mean(obj.spikes((condIds == i), :, :), 1));
            end
            obj.history{end+1} = ...
                sprintf('condition averaged. Number of condition : %d', n_cond);
        end
    end
end