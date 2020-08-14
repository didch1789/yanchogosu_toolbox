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

k = 1;
for i = 1:numel(divided)
    singlebatch = floor(total_length / numbatch);
    divided{i} = (k:singlebatch*i);
    
    if i == numel(divided)
        divided{i} = k:total_length;
    end
    
    k = k + singlebatch;
end


end
