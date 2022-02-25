function out = rmdupidx(input)
% input should be 2 column row, col indeces and function will remove same indeces.
% e.g)
%   input = [1 2; 2 1; 3 4]
%   out = rmdupidx(input)
%   out = [1 2; 3 4];
% 

if size(input, 2) ~=2 
    error('Input should be array with size #n X 2')
end

check1 = zeros(size(input));

for i = 1:size(input, 1)
    if input(i, 1) < input(i, 2)
        check1(i, :) = sort(input(i, :));
    end
end

out = unique(check1, 'rows');

end