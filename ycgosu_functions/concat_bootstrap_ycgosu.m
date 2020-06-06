function concat_bootstrap_ycgosu(datdirs, outputdir)
% :::Input:::
%   -datdirs: file directory of boots files (ex.{n X 1 filenames})
%   -outputdir file directory of output filename.
% :::Output:::
%   -bootout: WTS that contains conctatenated boot info.

w = [];

for i = 1:numel(datdirs)
    A = load(datdirs{i});
    fieldA = fieldnames(A);
    if numel(fieldA) > 1
        error('Your file should contain single field!!')
        break
    end
    w = [w;A.(fieldA{:}).w];
end

WTS = struct;

WTS.wste = squeeze(nanstd(w)); 
WTS.wmean = squeeze(nanmean(w)); 
WTS.wste(WTS.wste == 0) = Inf;  
WTS.wZ = WTS.wmean ./ WTS.wste;
WTS.wP = 2 * (1 - normcdf(abs(WTS.wZ)));

save(outputdir, 'WTS')
   

end