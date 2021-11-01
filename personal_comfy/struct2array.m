function out = struct2array(structs, dim)

if nargin < 2
    dim =1;
end

cellstructs = struct2cell(structs);
if dim == 1
    out = cat(1, cellstructs{:}); 
elseif dim == 2
    out = cat(2, cellstructs{:});
end
    

end