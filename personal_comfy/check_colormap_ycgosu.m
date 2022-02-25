function check_colormap_ycgosu(cmap, wide)
% cmap should be a n x 3 matrix
% wide should be a integer that indicates overall width of figure.(default:5)
if nargin < 2
    WidTh = 10;
else
    WidTh = wide;
end

x = zeros(50, size(cmap, 1)*WidTh,size(cmap, 2));
for i = 1:50
    x(i, :, :, :) = repelem(cmap, WidTh, 1);
end 
imshow(x)
set(gcf, 'Color', 'white')
end