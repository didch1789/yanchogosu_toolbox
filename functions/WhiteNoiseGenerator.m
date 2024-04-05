function [out, player] = WhiteNoiseGenerator(int_dur, int_freq, varargin)
% Inputs
%   int_sec:  How long do you want your output to be (in sec)
%   int_dur:  What frequency do you your output to be
% Outputs
%       out:  int_dur x int_freq @ 1 vector.
%    player:  Player object.  

doplay = false;
for i = 1:numel(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'magnitude', 'mag'}
                maginfo = varargin{i+1};
                % This should be 2 x n cell, where sum of the first row
                % should be the int_sec.
                % The second row should be the magnitudes.
                % e.g.) int_sec = 30, 
                %       maginfo = {[0 6], [6 12], [12 18], [18 24], [24 30],
                %                      1,  [1 2],    2   ,  [2 1],     1};
                %       % For mag, 1: default sound vol.
            case 'play'
                doplay     = true;
                playdevice = varargin{i+1};
                if isempty(playdevice)
                    playdevice = '';
                end
        end
    end
end

duration_is  = int_dur;  % in seconds
sample_rate  = int_freq;  % adjust as needed

% Calculate the number of samples
num_samples = duration_is * sample_rate;

% Generate white noise
white_noise = randn(1, num_samples);

% 
volume_ramp = ones(1, numel(white_noise));

if exist('maginfo', 'var')
    n_phase = size(maginfo, 2);
    for iP = 1:n_phase
        dur_time = maginfo{1, iP};
        ramps    = maginfo{2, iP};
        x_times  = (dur_time(1) * int_freq + 1):(dur_time(end) * int_freq);
        
        if isscalar(ramps)
            volume_ramp(x_times) = ramps;
        else
            volume_ramp(x_times) = linspace(ramps(1), ramps(end), numel(x_times));
        end

    end
end

increasing_volume_noise = white_noise .* volume_ramp;
out                     = increasing_volume_noise;

if doplay
    AudioInfo = audiodevinfo;
    AudioIdx  = contains({AudioInfo.output.Name}, playdevice);
    IDs       = cat(1, AudioInfo.output.ID);
    int_ID    = IDs(AudioIdx);
    
    player    = audioplayer(increasing_volume_noise, int_freq, 8, int_ID);
    play(player,[1 (get(player, 'SampleRate'))]);
end


end