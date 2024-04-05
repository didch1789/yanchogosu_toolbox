function out1 = clr_interpolate(in1, in2, nbin)

out1 = NaN(nbin, 3);
for ii = 1:3
    out1(:, ii) = linspace(in1(1, ii), in2(1, ii), nbin);
end

end