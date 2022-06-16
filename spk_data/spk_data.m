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
        
        function obj = normalize_spk(obj, varargin)
            rasterdats =  obj.spikes;
            for i = 1:numel(varargin)
                if ischar(varargin{i})
                    switch varargin{i}
                        case 'session_lvl'
                            norm_str = 'session';
                        case 'neuron_lvl'
                            norm_str = 'neuron';
                    end
                end
            end
            
            switch norm_str
                case {'session'}
                    mean_sess = cellfun(@(x) mean(x, 'all'), rasterdats, 'UniformOutput', false);
                    std_sess  = cellfun(@(x) std(x, 0, 'all'), rasterdats, 'UniformOutput', false);
                    obj.spikes = cellfun(@(x, y, z) (x-y)/z, rasterdats, mean_sess, std_sess, 'UniformOutput', false);
                case {'neuron'}
                    for iS = 1:numel(rasterdats)
                        [n_trial, n_neur, n_time] = size(rasterdats{iS});
                        singlesess       = permute(rasterdats{iS}, [2 3 1]);
                        singlesess_flat  = reshape(singlesess, n_neur, []);
                        singlesess_flat0 = zscore(singlesess_flat, 0, 2);
                        singlesess_0     = reshape(singlesess_flat0, n_neur, n_time, n_trial); 
                        rasterdats{iS}   = permute(singlesess_0, [3 1 2]);
                    end
            end
            obj.spikes = rasterdats;
            obj.history{end+1} = sprintf('normalized at %s lvl', norm_str);
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
        
        
        function obj = binning(obj, binsize, varargin)
            rasterdats = obj.spikes;
            for i = 1:numel(varargin)
                if ischar(varargin{i})
                    switch varargin{i}
                        case 'condavg'
                            rasterdats = obj.pseudopops.condavg;
                            obj.pseudopops.condavg = ...
                                cellfun(@(x) bin_raster_ycgosu(x, binsize), ...
                                rasterdats, 'UniformOutput', false);
                            obj.history{end+1} = ...
                                sprintf('binned on condavg with size %02d', binsize);
                            return;
                    end
                end
            end
            
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

        function obj = smoothing(obj, windowsize, varargin)
             rasterdats = obj.spikes;
            for i = 1:numel(varargin)
                if ischar(varargin{i})
                    switch varargin{i}
                        case 'condavg'
                            rasterdats = obj.pseudopops.condavg;
                            obj.pseudopops.condavg = ...
                                cellfun(@(x) smoothdata(x, 2, 'gaussian', windowsize), ...
                                rasterdats, 'UniformOutput', false);
                            obj.history{end+1} = ...
                                sprintf('smoothed on condavg with size %d', windowsize);
                            return;
                    end
                end
            end
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
        
        function out = return_condition_avg(obj)
            % would work after spk_data.condition_avg 
            n_cond = size(obj.pseudopops.condavg, 1);
            condavg_mat = cell2mat(obj.pseudopops.condavg')';
            n_time = size(condavg_mat, 1)/n_cond;
            out = mat2cell(condavg_mat, repmat(n_time, n_cond, 1), size(condavg_mat, 2));
        end
        
        
    end
end