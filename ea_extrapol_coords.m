function  ncoords_mm=ea_extrapol_coords(coords_mm,options)
% simple iterative extrapolation for more than 4 contacts
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn


ncoords_mm=zeros(options.elspec.numel,3);

for side=1
    ncoords_mm(1+(side-1)*options.elspec.numel:4+(side-1)*options.elspec.numel,:)=coords_mm(1+(side-1)*4:4+(side-1)*4,:);
for exc=5+(side-1)*options.elspec.numel:options.elspec.numel+(side-1)*options.elspec.numel
    ncoords_mm(exc,:)=ncoords_mm(exc-1,:)+(ncoords_mm(exc-1,:)-ncoords_mm(exc-2,:));
    
end

end

