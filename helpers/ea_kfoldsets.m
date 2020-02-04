function sets=ea_kfoldsets(n,k,random)
% outputs k discrete sets for CV applications looping over n data points

if ~exist('random','var')
    random=0;
end
if random
    allix=randperm(n);
else
    allix=1:n;
end

setlen=floor(n/k);
cnt=1;
for set=1:k
    if set==k
        sets{set}=allix(cnt:end);
    else
        sets{set}=allix(cnt:cnt+setlen-1);
    end
    cnt=cnt+setlen;
end




