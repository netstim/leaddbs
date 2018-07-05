function fitted_var=ea_fitto(var,sampledist)
% small helper function that will scale data to other distribution using
% a small GLM




if size(var,2)>size(var,1)
    var=var';
end
if size(sampledist,2)>size(sampledist,1)
    sampledist=sampledist';
end

b=glmfit(var,sampledist);
fitted_var=ea_addone(var)*b;



