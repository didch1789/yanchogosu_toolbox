function f = plot_raster_ycgosu(spikedat, varargin)
% INPUTS
%   spikedat: trial x timeseries of spike data of a single neuron (in ms)
%   varargin
%       'RasterColor' : rgb value. default color is black. (e.g. 'Color', [.5 .5 .5])
%       'AddLine'  : line index for indicating starting point. (e.g. 'Line',240)
%       'LineColor' : Color of the line. Also in RGB value. default color is red. (e.g. 'LineColor', [.5 .5 .5])
% OUTPUTS
%   f: figure handle

doline = false;
rastercolor = 'black';
linecolor = 'red';

for args = 1:numel(varargin)
    if isa(varargin{args}, 'char') || isa(varargin{args}, 'string')
        switch varargin{args}
            case 'RasterColor'
                rastercolor = varargin{args+1};
            case 'AddLine'
                doline = true;
                lineidx = varargin{args+1};
            case 'LineColor'
                linecolor = varargin{args+1};
            otherwise
                continue
        end
        
    end
end


%% Raster
spikes = spikedat;
trials = size(spikes, 1);
triallength = size(spikes, 2);

for i = 1:trials
    spks = find(spikes(i, :));
    xspks = repmat(spks, 2, 1);
    yspks = nan(size(xspks));
    
    yspks(1, :) = i - 1;
    yspks(2, :) = i;
    
    plot(xspks, yspks, 'Color', rastercolor);hold on
end
set(gca, 'XLim', [0, triallength], 'Ylim', [0, trials], 'TickDir', 'out', 'TickLength', [.005 .005])

if doline
   line([lineidx lineidx], [0 trials], 'Color', linecolor)
end

hold off;

f = gcf;

end
