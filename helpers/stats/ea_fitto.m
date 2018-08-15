function [fitted_var,b,stats,Rsquare,Fstat,poverall]=ea_fitto(var,sampledist)
% small helper function that will scale data to other distribution using
% a small GLM




if size(var,2)>size(var,1)
    var=var';
end
if size(sampledist,2)>size(sampledist,1)
    sampledist=sampledist';
end

[b,dev,stats]=glmfit(var,sampledist);
fitted_var=ea_addone(var)*b;

[~,~,~,~,statr]=regress(sampledist,ea_addone(var));

Rsquare=statr(1);
Fstat=statr(2);
poverall=statr(3);

