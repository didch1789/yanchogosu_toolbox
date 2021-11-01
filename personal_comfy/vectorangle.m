function degreeresults = vectorangle(u, v)
% x and y should be both n x 1 vector.
% output yields only 0 - 360 degree

CosTheta = dot(u,v) / (norm(u) * norm(v));
degreeresults = real(acosd(CosTheta));


end

