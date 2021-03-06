function circle_plot_ycgosu(names, values, varargin)
% plot your characters in different size in circle wise manner according to their assigned number. 
% inputs:
%       names: cell array of names to be plotted.
%       values: same size with names. should be a array with numbers.
% varagin: % not yet...!
%       'dosize': will draw 
%       'docolor':
%       'BackColor':
%       'Colormap': 
% 

dosize = 1;
docolor = 0;
backcolor = [1 1 1];
scalefactor = numel(names);
posmap = [254,224,210;252,187,161;252,146,114;251,106,74;239,59,44;203,24,29;165,15,21];
negmap = [222,235,247;198,219,239;158,202,225;107,174,214;66,146,198;33,113,181;8,81,156];

len = size(names, 1);
x = linspace(0, 360-(360/len), len);

[~, sortidx_temp] = sort(values);
[~, sortidx] = sort(sortidx_temp);
scales = linspace(500/scalefactor, 1200/scalefactor, len);
scales = scales(sortidx);
fontsize = 500/scalefactor;

for i=1:numel(varargin)
    if ischar(varargin{i})
        switch varargin{i} 
            case {'dosize'}
                dosize = 1;     docolor = 0;
            case {'docolor'}
                docolor = 1;    dosize = 0;
            case {'BackColor'}
                backcolor = varargin{i+1};
            case {'posnegmap'}
                posmap = varargin{i+1};
                negmap = varargin{i+2};
            case {'fontscale'}
                scalefactor = varargin{i+1};
            otherwise
                error('Unidentified varargin') 
        end
    end
end

if dosize
    
    for i = 1:len
        scatter(cosd(x(i)), sind(x(i)),  'MarkerFaceColor', backcolor, 'MarkerEdgeColor', backcolor)
        if (x(i) < 90 || x(i) > 270)
            text(cosd(x(i)), sind(x(i)), names(i), ...
                'FontSize', scales(i), 'Rotation', x(i), 'Interpreter', 'none');hold on;
        else
            text(cosd(x(i)), sind(x(i)), names(i), 'FontSize', scales(i),...
                'Rotation', x(i)+180, 'HorizontalAlignment', 'right', 'Interpreter', 'none');hold on;
        end
    end
    xlim([-2 2]);ylim([-2 2])
    set(gcf, 'Color', backcolor)
    set(gca, 'Color', backcolor, 'Xcolor', 'none', 'Ycolor', 'none')
    
elseif docolor
    [sortedvals, sortidx_temp] = sort(values, 'descend');
    sorted_field = names(sortidx_temp);
    numpos = sum(sortedvals > 0);
    numneg = sum(sortedvals < 0);

    redmap = zeros(numpos, 3);bluemap = zeros(numneg, 3);

    for col = 1:3
        redmap(:, col) = linspace(posmap(1, col), posmap(end, col), numpos);
        bluemap(:, col) = linspace(posmap(1, col), negmap(end, col), numneg);
    end
    
    redmap = redmap(end:-1:1, :) ./ 256;
    bluemap = bluemap ./ 256;
        
    R = 1; B = 1;
    for i = 1:len
        scatter(cosd(x(i)), sind(x(i)), 'MarkerFaceColor', backcolor, 'MarkerEdgeColor', backcolor)
        if sortedvals(i) > 0 && (x(i) < 90 || x(i) > 270)
            text(cosd(x(i)), sind(x(i)), sorted_field(i), 'Color', redmap(R, :), ...
                'FontSize', fontsize, 'Rotation', x(i), 'Interpreter', 'none');hold on;
            R = R + 1;
        elseif sortedvals(i) > 0 && ~(x(i) < 90 || x(i) > 270)
            text(cosd(x(i)), sind(x(i)), sorted_field(i), 'Color', redmap(R, :), ...
                'FontSize', fontsize, 'Rotation', x(i) + 180, 'Interpreter', 'none', 'HorizontalAlignment', 'right');hold on;
            R = R + 1;
        elseif sortedvals(i) < 0 && (x(i) < 90 || x(i) > 270)
            text(cosd(x(i)), sind(x(i)), sorted_field(i), 'Color', bluemap(B, :), ...
                'FontSize', fontsize, 'Rotation', x(i), 'Interpreter', 'none');hold on;
            B = B + 1;
        elseif sortedvals(i) < 0 && ~(x(i) < 90 || x(i) > 270)
            text(cosd(x(i)), sind(x(i)), sorted_field(i), 'Color', bluemap(B, :), ...
                'FontSize', fontsize, 'Rotation', x(i)+180, 'Interpreter', 'none', 'HorizontalAlignment', 'right');hold on;
            B = B + 1;
        end
    end
    xlim([-2 2]);ylim([-2 2])
    set(gcf, 'Color', backcolor)
    set(gca, 'Color', backcolor, 'Xcolor', 'none', 'Ycolor', 'none')
end



end