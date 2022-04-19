classdef spk_data
    properties
        spikes = {};
        pseudopops = struct;
        dims = 'trial X neuron X time'
        history = {'Processed History:'};
    end
    
    methods
        
        function obj = spk_more_than(obj, spkthr, varargin)
            % spkthr      : unit is spk/sec.
            % interval : unit is ms.
            spkinterval = 1;
            for vi = 1:numel(varargin)
                if ischar(varargin{vi})
                    switch varargin{vi}
                        case 'interval'
                            spkinterval = varargin{i+1};
                    end
                end
            end
            rasterdats =  obj.spikes;
            for i = 1:numel(rasterdats)
                rasterdat = rasterdats{i};
                if numel(size(rasterdat)) > 2
                    avgraster = mean(rasterdat, [1, 3]);
                    time_range = size(rasterdat, 3);
                else
                    avgraster = mean(rasterdat, [1, 2]);
                    time_range = size(rasterdat, 2);
                end
                spk_per_sec = avgraster ./ (time_range/spkinterval) * 1000;
                neuron_idx = spk_per_sec > spkthr;
                if numel(size(rasterdat)) > 2
                    obj.spikes{i} = rasterdat(:, neuron_idx, :);
                else
                    if neuron_idx == true
                        obj.spikes{i} = rasterdat;
                    else
                        obj.spikes{i} = NaN;
                    end 
                end
            end
            obj.history{end+1} = sprintf('only spikes more than %d per/sec', spkthr);
        end
        
        function obj = zscore_spk(obj, dim)
            rasterdats =  obj.spikes;
            obj.spikes = cellfun(@(x) zscore(x, 0, dim), rasterdats, 'UniformOutput', false);
            obj.history{end+1} = sprintf('zscored in %d dim', dim);
        end
       
        
        function obj = time_range(obj, t_range)
            rasterdats =  obj.spikes;
            for i = 1:numel(rasterdats)
                rasterdat = rasterdats{i};
                if numel(size(rasterdat)) > 2 
                    obj.spikes{i} = rasterdat(:, :, t_range);
                else
                    obj.spikes{i} = rasterdat(:, t_range);
                end
            end
            obj.history{end+1} = sprintf('int_time:%d-%d', t_range(1), t_range(end));
        end
        
        
        function obj = binning(obj, binsize)
            rasterdats = obj.spikes;
            for i = 1:numel(rasterdats)
                rasterdat = rasterdats{i};
                if numel(size(rasterdat)) > 2 
                    n_trial  = size(rasterdat, 1);
                    n_neuron = size(rasterdat, 2);
                    n_time   = size(rasterdat, 3);
                    Q        = fix(n_time / binsize);
                    R        = mod(n_time, binsize);
                    binned   = NaN(n_trial, n_neuron, Q);

                    for tt = 1:n_trial
                        spktrials   = rasterdat(tt, :, :);
                        l = 1;
                        for bb = 1:Q
                            if R ~= 0
                                k = binsize + 1;
                                R = R - 1;
                            else
                                k = binsize;
                            end
                            binned(tt, :, bb) = squeeze(sum(spktrials(:, :, l:l+k-1), 3));
                            l = l + k;
                        end
                    end
                else
                    binned = bin_raster_ycgosu(rasterdat, binsize);
                end
                obj.spikes{i} = binned;
            end
            obj.history{end+1} = sprintf('binned with size %02d', binsize);
        end

        function obj = smoothing(obj, windowsize)
            rasterdats = obj.spikes;
            for i = 1:numel(rasterdats)
                rasterdat = rasterdats{i};
                if numel(size(rasterdat)) > 2
                    obj.spikes{i} = smoothdata(rasterdat, 3, 'gaussian', windowsize);
                else
                    obj.spikes{i} = smoothdata(rasterdat, 2, 'gaussian', windowsize);
                end
            end
            obj.history{end+1} = sprintf('smoothed with size %d', windowsize);
        end
        
        

        function obj = condition_avg(obj, condIds, varargin)
            % "condIds" should be cell.
            for i = 1:numel(varargin)
                if ischar(varargin{i})
                    switch varargin{i}
                        case 'condnames'
                            obj.pseudopops.condnames = varargin{i+1};
                        case 'n_sample'
                            n_sample = varargin{i+1};
                    end
                end
            end
            
            n_sess = numel(obj.spikes);
            
            for sidx = 1:n_sess
                d_spk  = obj.spikes{sidx};      d_cond = condIds{sidx};
                saveidx = ~isnan(d_cond);
                
                d_cond = d_cond(saveidx);       d_spk  = d_spk(saveidx, :, :);
                n_cond = numel(unique(d_cond));
                
                for cidx = 1:n_cond
                    int_idx = find(d_cond == cidx);
                    if exist('n_sample', 'var')
                        sample_idx = datasample(int_idx, n_sample, 'Replace', false);
                    else
                        sample_idx = int_idx;
                    end
                    c_spk = mean(d_spk(sample_idx, :, :), 1);
                    obj.pseudopops.condavg{cidx, sidx} = permute(c_spk, [2 3 1]);
                    obj.pseudopops.trialnum(cidx, sidx)= sum(d_cond == cidx);
                end

            end
            
            obj.history{end+1} = ...
                sprintf('condition averaged. Number of condition : %d', n_cond);
        end
        
        
    end
end