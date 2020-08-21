function orthviews_fsl(obj, varargin)
% Similar to orthviews but shows brain map in fsleyes.
% you should have fsl in you "/usr/local/" (symbolic link will do.)
% Makes temporary .nii file in current directoty since fsl cannot read fmri_data or stat_img
% directly. It will be deleted soon after fsleyes shows up.

% input:
%   obj: fmri_data or statisting image object. 
%   varargin:
%       'underlay'                 : directory of underlay img.
%                                    (default: keuken_enhanced_for                            
%       'threshold_val', [num num] : only display values in num~num;
%       'threshold_p', [num]       : only display values below p.
%                                    (only works in statsting image) 
%       'colormap', string         : refer to fsleyes colormap(default: 'hot')
%       'auto_thresh'              : will display top and bottom 95 - 99 %
%                                    values, 
%       'only_pos'                 : will display only positive values
%       'only_neg'                 : will display only negative values
cmap = 'hot';
underlay = which('keuken_2014_enhanced_for_underlay.nii');

for i = 1:numel(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'underlay'
                underlay = varargin{i+1};
            case 'threshold_val'
                threshold_val = varargin{i+1};
            case 'threshold_p'
                threshold_p = varargin{i+1};
            case 'colormap'
                cmap = varargin{i+1};
            case 'auto_thresh'
                auto_thresh = true;
            case 'only_pos'
                only_pos = true;
            case 'only_neg'
                only_neg = true;
        end
    end
end

obj.fullpath = fullfile(pwd, 'for_orthivews_fsl.nii');
obj.write;
dat = obj.dat;
b1 = prctile(dat, 1);   b2 = prctile(dat, 5);
a1 = prctile(dat, 95);  a2 = prctile(dat, 95);


if 
    cmd = sprintf('fsleyes %s %s -dr %.2f %.2f -cm %s',  cmap)


end