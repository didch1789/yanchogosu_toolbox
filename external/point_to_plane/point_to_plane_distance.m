function [d2H, H] = point_to_plane_distance(M, n, I)
%% point_to_plane_distance : function to compute the distance between the 3D point M
% and the plane (I,n). Also provides the coordinates of H,
% the projection of M on (I,n), and also works for a list of points.
%
% Author and support nicolas.douillet (at) free.fr, 2019-2020.
%             
%
% Syntax
%
% d2H = point_to_plane_distance(M, n, I);
% [d2H, H] = point_to_plane_distance(M, n ,I);
%
%
% Description
%
% d2H = point_to_plane_distance(M, n, I) computes the distance d2H between M and the
% plane (I,n).
%
% [d2H, H] = point_to_plane_distance(n ,I, M) also provides the coordinates
% of H, the projection of M on (I,n).
%
%
% Input arguments
%
%       [ |  |  |]
% - M = [Mx My Mz], real row vector -or matrix- double , the point -or list of N points- which
%       [ |  |  |]
%
%       we want the distance to the plane. size(M) = [N,3].
%
% - n = [nx ny nz], real row vector double, a plane normal vector, size(n) = [1,3].
%
% - I = [Ix Iy Iz], real row vector double, a point belonging to the plane, size(I) = [1,3]. 
% 
%
% Output arguments
%
% - d2H : real scalar -or vector- double, the euclidian distance between M and the plane (I,n). size(d2H) = [N,1].
%
% - H : real vector -or matrix- double, the projected point(s) on the plane. size(H) = [N,3].
%
%
% Example #1
% Single point
%
% M = [1 1 1];
% I = [1 0 0];
% n = M;
% [d2H, H] = point_to_plane_distance(M, n ,I) % expected distance : 2/sqrt(3); expected coordinates : [1/3 1/3 1/3]                                                                                                            
%
%
% Example #2
% List of points
% M = [1 1 1; 2 1 1];
% I = [1 0 0];
% n = [1 1 1];
% Expected distances : 2/sqrt(3), sqrt(3); expected coordinates : [1/3 1/3 1/3], [1 0 0].
% [d2H, H] = point_to_plane_distance(M, n ,I)


%% Inputs parsing
assert(nargin > 2, 'Not enough input arguments.');
assert(nargin < 4, 'Too many input arguments.')

nb_pts = size(M,1);

assert(isequal(size(n),size(I),[1,3]) || ...        
       isequal(ndims(M),ndims(n),2),...
       'Inputs arguments n and I must have the same format (size, number of elments, and be of dimension 2).');

assert(isequal(size(I,2),size(M,2),3),'Input arguments M, n, and I must have the same number of columns (3).');   
assert(isreal(M) && isreal(n) && isreal(I),'Input arguments M, n, and I must be real.');


%% Body
d_I = -(n(1)*I(1)+n(2)*I(2)+n(3)*I(3));
t_H = -(repmat(d_I,[nb_pts,1])+n(1)*M(:,1)+n(2)*M(:,2)+n(3)*M(:,3)) / sum(n.^2);

x_H = M(:,1) + t_H*n(1);
y_H = M(:,2) + t_H*n(2);
z_H = M(:,3) + t_H*n(3);

H = cat(2,x_H,y_H,z_H);

d2H = sqrt(sum((M-H).^2,2));


end % point_to_plane_distance