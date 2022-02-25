function statimg = make_statimg_ycgosu(fmriobj, WTS)
% Inputs:
%       fmriobj: stats.weight_obj of fmri_data.predict
%       WTS: struct that contains wste, wmean, wZ, wP
% Output:
%       statimg = statimg without contrast

M = statistic_image;
M.p = WTS.wP';
M.ste = WTS.wste';
M.dat = fmriobj.dat;
M.volInfo = fmriobj.volInfo;

statimg = M;

end