function output = xforxin(fromm, till, str)
% inspired by [x for x in range(10)] in Python!
% output gives cell array of 'str fromm' ~ 'str till'
% e.g)
% xforxin(2, 10, 'num')
% {'run02' }
w = floor(log10(max(fromm, till)));
output = cell(abs(till-fromm+1), 1);

k = 1;
for i = fromm:till
    output{k} = sprintf([str, '%0' num2str(w+1), 'd'], i);
    k = k + 1;
end

end