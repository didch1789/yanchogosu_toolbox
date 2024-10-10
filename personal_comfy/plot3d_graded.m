function out_axis = plot3d_graded(XYZ, clrs, varargin)

linewidth = 1;
linespec = '-';
for i = 1:numel(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'LineWidth'}
                linewidth = varargin{i+1};
            case {'LineSpec'}
                linespec = varargin{i+1};
        end
    end

end

n_time = size(XYZ, 1);
n_clr  = size(clrs, 1);

if ~(n_time == n_clr)
    error('number of color input does not match with that of the XYZ input')
end

for t_i = 1:n_time
    if t_i == n_time
        break
    end
    tidx = t_i:(t_i+1); 

    x = XYZ(tidx, 1);
    y = XYZ(tidx, 2);
    z = XYZ(tidx, 3);

    plot3(x, y, z, 'Color', clrs(t_i, :), 'LineWidth', linewidth, 'LineStyle', linespec); hold on;
end

out_axis = gca;


end