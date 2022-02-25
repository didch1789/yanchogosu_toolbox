function out = make_cvparition(n, varargin)

kfold = 5;
for i = 1:numel(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'kfolds'
                kfold = varargin{i+1};
            otherwise
                error('No such option')
        end
    end
end

random_n = randperm(n);
sets = divide_ycgosu(n, kfold);

for i = 1:numel(sets)
    out.teIdx{i} = random_n(sets{i});
    out.trIdx{i} = random_n(~ismember(random_n, random_n(sets{i})));
end
            
                
        

end