function check_colormap_ycgosu(cmap)
% cmap should be a n x 3 matrix
x = zeros(10, size(cmap, 1),size(cmap, 2));
for i = 1:10
    x(i, :, :, :) = cmap;
end 
imshow(x)
end