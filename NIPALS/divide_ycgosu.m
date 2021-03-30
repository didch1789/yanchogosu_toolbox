function divided = divide_ycgosu(total_length, numbatch)

%%% inputs
% total_length: total number(can be a double or whole numbers
% numbatch: number of batches you want to divide

%%% output
% divided number in cell type

if numel(total_length) ~= 1
    total_length = max(total_length);
end

divided = cell(numbatch, 1);
Q = fix(total_length/numbatch);
R = rem(total_length,numbatch);
k = 0;
for i = 1:numel(divided)
    if R ~= 0
        divided{i} = (Q*(i-1)+1+k:Q*i+1+k);
        R = R - 1;
        k = k + 1;
    else
        divided{i} = (Q*i-(Q-1):Q*i);
    end
end


end
