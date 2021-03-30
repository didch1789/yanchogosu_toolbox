function [out desiign] = ActResttime(wholeTRs, actonsets, varargin)
% % INPUTS % 
%
% wholeTRs  : number of intTR (e.g. 240)
% actonsets : stimulus onset in seconds (e.g. [5 9 10 15])
% TR(in sec): TR in your study (e.g. 2)
%
% varargin (weighted mean or un-weighted mean)
%
% or you are to use un-weighted mean, specify "delay" and "windows"
%   durs  : duration of the stimulus (in sec)
%   windows   : How many TRs would you want get from the delay? 
%               (e.g. [2 2] would mean that you are extracting 2 TRs before and after the peak of HRF)
%
% (e.g. ActResttime(180, [6 20 40 50 70 100 200], 'TR', 2, 'duration', 4, 'windows', [2 3])
%
% % OUTPUTS %
% out.actTR  : string of TR of activation
% out.restTR : string of TR of rest. (whole TR - act)

whatchar = ':';

for v = 1:numel(varargin)
    if isa(varargin{v}, 'char')
        switch varargin{v}
            case 'TR'
                TR = varargin{v+1};
            case  'duration'
                durs = varargin{v+1};            
            case  'windows'
                windows = varargin{v+1};            
            case 'joints'
                whatchar = varargin{v+1};          
        end
    end
end

OnsDurinput = [actonsets repmat(durs, size(actonsets, 1), 1)];
HRFfunc = onsets2fmridesign(OnsDurinput, TR, TR*wholeTRs, spm_hrf(TR));

[~, peakidx] = findpeaks(HRFfunc(:, 1));

desiign = HRFfunc(:, 1);

actTR = '';
if whatchar == '-'
    restTR = '0';
else
    restTR = '1';
end

for i = 1:numel(peakidx) 
    st = peakidx(i)-windows(1);
    ed = peakidx(i)+windows(2);
    if whatchar == '-'
        st = st-1;
        ed = ed-1;
    end
    restTR = [restTR, whatchar ,num2str(st-1), ','];
    actTR = [actTR, num2str(st), whatchar, num2str(ed), ','];
    restTR = [restTR, num2str(ed+1)];
    if i == numel(peakidx) 
        actTR(end) = [];
        if whatchar == '-'
            restTR = [restTR, whatchar, num2str(wholeTRs-1)];
        else
            restTR = [restTR, whatchar, num2str(wholeTRs)];
        end
    end
end



out.actTR = actTR;
out.restTR = restTR;







end