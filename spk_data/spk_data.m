classdef spk_data
    properties
        spikes = {};
        pseudopops = struct;
        dims = 'trial X neuron X time'
        history = {'Processed History:'};
    end
    
    methods
        
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
                    end
                end
            end
            
            n_sess = numel(obj.spikes);
            
            for sidx = 1:n_sess
                d_spk  = obj.spikes{sidx};      d_cond = condIds{sidx};
                n_cond = numel(unique(d_cond));
                
                for cidx = 1:n_cond
                    obj.pseudopops.condavg{sidx, cidx} = squeeze(mean(d_spk(d_cond == cidx, :, :), 1));
                    obj.pseudopops.trialnum(sidx, cidx)= sum(d_cond == cidx);
                end

            end
            
            obj.history{end+1} = ...
                sprintf('condition averaged. Number of condition : %d', n_cond);
        end
        
        
    end
end