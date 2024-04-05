function divided = divide_ycgosu(total_length, numcell)

%%% inputs
% total_length: total number(can be a double or whole numbers
% numbatch: number of batches you want to divide

%%% output
% divided number in cell type

if numel(total_length) ~= 1
    total_length = max(total_length);
end


divided = cell(numcell, 1);
Q = fix(total_length/numcell);
R = rem(total_length,numcell);
k = 0;
for i = 1:numel(divided)
    if R ~= 0
        divided{i} = (Q*(i-1)+k+1:Q*i+k+1);
        R = R - 1;
        k = k + 1;
    else
        divided{i} = (Q*(i-1)+k+1:Q*(i)+k);
    end
end  
    


end
