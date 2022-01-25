%% point_to_plane_distance
%
% Function to compute the distance between the 3D point M
% and the plane (I,n). Also provides the coordinates of H,
% the projection of M on (I,n), and also works for a list of points.
%
% Author and support nicolas.douillet (at) free.fr, 2019-2020.
%                                                                       
%% Syntax
%
% d2H = point_to_plane_distance(M, n, I);
%
% [d2H, H] = point_to_plane_distance(M, n ,I);
%
%% Description
%
% d2H = point_to_plane_distance(M, n, I) computes the distance d2H between M and the
% plane (I,n).
%
% [d2H, H] = point_to_plane_distance(M, n ,I) also provides the coordinates
% of H, the projection of M on (I,n).
%
%% See also
%
% <https://fr.mathworks.com/matlabcentral/fileexchange/73478-point-to-line-distance-3d-2d?s_tid=prof_contriblnk point_to_line_distance> |
% <https://fr.mathworks.com/matlabcentral/fileexchange/73760-line-plane-intersection-3d?s_tid=prof_contriblnk line_plane_intersection> |
% <https://fr.mathworks.com/matlabcentral/fileexchange/73853-planes-intersection planes_intersection> |
% <https://fr.mathworks.com/matlabcentral/fileexchange/73922-lines-intersection-3d-2d?s_tid=prof_contriblnk lines_intersection>
%
%% Input arguments
%
%        [ |  |  |]
% - M = [Mx My Mz], real row vector -or matrix- double , the point -or list of N points- which
%        [ |  |  |]
%
%       we want the distance to the plane. size(M) = [N,3].
%
% - n = [nx ny nz], real row vector double, a plane normal vector, size(n) = [1,3].
%
% - I = [Ix Iy Iz], real row vector double, a point belonging to the plane, size(I) = [1,3]. 
% 
%% Output arguments
%
% - d2H : real scalar -or vector- double, the euclidian distance between M and the plane (I,n). size(d2H) = [N,1].
%
% - H : real vector -or matrix- double, the projected point(s) on the plane. size(H) = [N,3].

%% Example #1
% Single point

M = [1 1 1];
I = [1 0 0];
n = M;
[d2H, H] = point_to_plane_distance(M, n ,I) % expected distance : 2/sqrt(3); expected coordinates : [1/3 1/3 1/3]                                                                                                            

%% Example #2
% List of points

M = [1 1 1; 2 1 1];
I = [1 0 0];
n = [1 1 1];
% Expected distances : 2/sqrt(3), sqrt(3); expected coordinates : [1/3 1/3 1/3], [1 0 0].
[d2H, H] = point_to_plane_distance(M, n ,I)