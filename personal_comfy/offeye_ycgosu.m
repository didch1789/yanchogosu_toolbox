function out = offeye_ycgosu(num)
    % put 1s in inversed diagonal space.
    % e.g) offeye_ycgosu(3)
%             =  [0 0 1
%                 0 1 0
%                 1 0 0];
    if size(num, 1) == size(num, 2) || num(1) == num(2)
        sz = num;
        out = eye(sz);
        out = out(:, end:-1:1);
    else
        row = num(1) ; col = num(2);
        if row > col
            out = eye(col);
            out = out(:, end:-1:1);
            out = [out zeros(size(out, 1), row-col)];
        elseif row < col
            out = eye(row);
            out = out(:, end:-1:1);
            out = [out; (zeros(size(out, 2), col - row))'];
        end
    end
    
    
end