function nsols = numevalsol(sols,varargin)

for k = 1:(nargin-1)/2,
    eval(sprintf('%s = varargin{%i}; ',varargin{2*k-1},2*k));
end;
nsols = struct;
fn = fieldnames(sols);
 strr = @(x) strrep(strrep(strrep(x,'*','.*'),'^','.^'),'/','./');

for k = 1:length(fn),    
    f = getfield(sols,fn{k});
    clear res
    for j = 1:length(f);
        res(j,:) = eval(strr(char(f(j))));
    end;
    nsols = setfield(nsols,fn{k},res);
    
end;