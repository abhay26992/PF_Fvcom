% Merging structure fields having the same field names
% Written by Abhay Prakash (abhay.prakash@natgeo.su.se)

function S = CatStructFields(varargin)
F = cellfun(@fieldnames,varargin,'uni',0);
assert(isequal(F{:}),'All structures must have the same field names.')
T = [varargin{:}];
S = struct();
F = F{1};
for k = 1:numel(F)
    tv=T.(F{k}); 
    l_tv=length(size(tv));
    if l_tv==3 %if number of dimensions = 3 (time = 3rd dim)
        dim=3;
        S.(F{k}) = cat(dim,T.(F{k}));
    elseif l_tv==2 %if number of dimensions = 2 (time = 2nd dim)
        dim=2;
        S.(F{k}) = cat(dim,T.(F{k}));
    end
end
end