function outangles = angle_distrib_between_planes(plane1, plane2, varargin)
% calculate random 5000 (default) vectors in plane1 column space and plane2 column space.
% Angle value is in between 0 ~ 90 degree.
% varargin
%   'numperm' : number of random permutation. 
% seems no big use...
numperm = 5000;
for vi = 1:numel(varargin)
    if ischar(varargin{vi})
        switch varargin{vi}
            case 'numperm'
                numperm = varargin{vi+1};
        end
    end
end
[ndim1, nvec1] = size(plane1);
[ndim2, nvec2] = size(plane2);

if ndim1 ~= ndim2
    error('Check dimensions of the input!')
end

outangles = NaN(numperm, 1);
for i = 1:numperm
    weight1 = orth(randn(nvec1, 1));
    weight2 = orth(randn(nvec2, 1));

    colsp1 = (plane1 * weight1);
    colsp2 = (plane2 * weight2);
    
    outangles(i) = corr(colsp1, colsp2); 
end


end