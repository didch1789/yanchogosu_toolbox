function extract_statistics(inputobj, varargin)
% make beta images in the current directory.

for i = 1:numel(varargin)
    if isa(varargin{i}, 'char')
        switch varargin{i}
            case 'prefix'
                prefix = varargin{i+1};
        end
    end
end

prfx = '';

if exist('prefix', 'var'), prfx = [prefix, '_'];end

template = inputobj;
type = [prfx,template.type];
outdir = sprintf('%smaps_yc', type);
if ~exist(outdir, 'dir'), mkdir(outdir);end

for i = 1:size(inputobj.p, 2)
    template.p          = inputobj.p(:, i);
    template.ste        = inputobj.ste(:, i);
    template.threshold  = inputobj.sig(:, i);
    template.dat        = inputobj.dat(:, i);
    template.removed_images = inputobj.removed_images(i);
    template.fullpath = fullfile(pwd, outdir, sprintf('%s%04d.nii', template.type, i));
    write(template, 'overwrite');
end

    
    
