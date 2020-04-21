function divided = divide_ycgosu(total_length ,interval)

%%% inputs
% total_length: total number(can be a double or whole numbers
% interval: numbers to be divided in one interval

%%% output
% divided number in cell type

if numel(total_length) ~= 1
    total_length = max(total_length);
end

divided = {};
k = 1;
l = 1;
for i = 1:total_length
    
    if mod(i, interval) == 0
        divided{l} = [k:i];
        l = l + 1;
        k = k + interval;
    end
    
    if mod(total_length, interval) ~= 0
        divided{l} = [k:k+mod(total_length, interval) - 1];
    end
    
end

end
