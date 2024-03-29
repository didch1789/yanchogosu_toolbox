function [f, b, plotdata] = plot_psth_ycgosu(spikedat, nbins, varargin)
% INPUTS
%
%   spikedat: trial x timeseries of spike data of a single neuron (in ms)
%
%   nbins: how many bins you want in your PSTH. (e.g. 'nbins', 10)
%          If your input time series are 1000ms long and define binsize as 100, each bin will contain 10 ms timeseries data. 
%            
%   varargin
%       'FaceColor' : FaceColor of histogram. default color is grey (ref. histogram)
%       'FaceAlpha' : FaceAlpha of histogram. default alpha is 0.8 (ref. histogram)
%       'EdgeColor' : EdgeColor of histogram. default color is white (ref. histogram)
%       'AddLine'   : Draw vertical line at defined index (e.g. 'AddLine', 250)
%       'LineColor' : Color of the vertical line. default color is red.
%
% OUTPUTS
%   f: figure handle
%   b: bar handle
%   plotdata : data plotted on the PSTH.
%
% @Jungwoo. 2021.07.22

if nargin < 2
   error('Not enough inputs!!') 
end


facecolor   = [.3 .3 .3];
facealpha   = .8;
edgecolor   = 'white';
doline      = false;
linecolor   = 'red';


for args = 1:numel(varargin)
    if isa(varargin{args}, 'char') || isa(varargin{args}, 'string')
        switch varargin{args}
            case 'FaceColor'
                facecolor = varargin{args+1};
            case 'FaceAlpha'
                facealpha = varargin{args+1};
            case 'EdgeColor'
                edgecolor = varargin{args+1};
            case 'AddLine'
                doline = true;
            case 'LineColor'
                linecolor = varargin{args+1};
            otherwise
                continue
        end
        
    end
end


%% PSTH
spikes = spikedat;
trials = size(spikes, 1);
triallength = size(spikes, 2);

allspks = [];
for i = 1:trials
    spks = find(spikes(i, :));
    allspks = [allspks;spks'];
end

bindur = triallength / nbins;
nobins = 1000 / bindur;                         % number of bins/sec

hcounts = histcounts(allspks, nbins);
unihcounts = unique(hcounts);

for yt = 1:numel(unihcounts)
   oldy = unihcounts(yt);
   conv = (oldy / trials) * nobins;             % convert to Hz
   hcounts(hcounts == oldy) = round(conv, 2);
end

plotdata = hcounts;

b = bar((0+(bindur)/2):bindur:bindur*nbins, hcounts);
set(b, 'FaceColor', facecolor, 'FaceAlpha', facealpha, 'EdgeColor', edgecolor)
set(gca, 'XLim', [0, triallength],'TickDir', 'out', 'TickLength', [.005 .005])

if doline
   line([lineidx lineidx], [0 max(h.Values)], 'Color', linecolor);
end

f = gcf;

end
