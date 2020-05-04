function out = bootstrp_ycgosu(dats, saveornot)
% dats : construct that contains bootstrapped weights in field name 'w' 
% saveornot[1 0]: 1 = save boot weights, 0 = delete boot weights
WTS.w = dats.w;

% got it from fmri_data.predict
WTS.wste = squeeze(nanstd(WTS.w)); %1/20/16 add squeeze for multiclass case
WTS.wmean = squeeze(nanmean(WTS.w)); %1/20/16 add squeeze for  multiclass case
WTS.wste(WTS.wste == 0) = Inf;  % in case unstable regression returns all zeros
WTS.wZ = WTS.wmean ./ WTS.wste;  % tor changed from wmean; otherwise bootstrap variance in mean inc in error; Luke renamed to avoid confusion
WTS.wP = 2 * (1 - normcdf(abs(WTS.wZ)));

if saveornot == 0
    WTS.w = [];
end 
    
out = WTS;

end