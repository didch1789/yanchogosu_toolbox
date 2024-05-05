function [lineobj, patchobj] = plotShadedErrorbar(x, y, varargin)

if size(x, 1) > size(x, 2), x = x';end
if size(y, 1) > size(y, 2), y = y';end

for i = 1:numel(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'LineColor'}
                lc  = varargin{i+1};
            case {'LineWidth'}
                lw  = varargin{i+1};
            case {'ShadeColor'}
                sc  = varargin{i+1};
            case {'ShadeEdgeColor'}
                sec = varargin{i+1};
            case {'ShadeAlpha'}
                sa  = varargin{i+1};
            case {'Error'}
                err = varargin{i+1};
                if size(err, 1) > size(err, 2), err = err';end
        end
    end
end


lineobj = plot(x, y, 'Color', lc, 'LineWidth', lw);hold on
patchobj = patch([x flip(x)], [y - err, flip(y+err)], sc, 'FaceAlpha',sa, 'EdgeColor', sec);



end