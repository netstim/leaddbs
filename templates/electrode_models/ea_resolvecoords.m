function [coords,trajectory]=ea_resolvecoords(markers,options)
% This function inflates the 2 markers reconstructed in
% ea_reconstruct_trajectory to the contacts as specified by the lead model.
% __________________________________________________________________________________
% Copyright (C) 2015 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

load([options.earoot,'templates',filesep,'electrode_models',filesep,options.elspec.matfname]);
for side=1:length(markers)
M=[markers(side).head,1;markers(side).tail,1;markers(side).x,1;markers(side).y,1];
E=[electrode.head_position,1;electrode.tail_position,1;electrode.x_position,1;electrode.y_position,1];
X=linsolve(E,M);
coords_mm=[electrode.coords_mm,ones(size(electrode.coords_mm,1),1)];
coords{side}=X'*coords_mm';
coords{side}=coords{side}(1:3,:)';






trajvector{side}=(markers(side).tail-markers(side).head)/norm(markers(side).tail-markers(side).head);
trajectory{side}=[markers(side).head-trajvector{side}*5;markers(side).head+trajvector{side}*25];
trajectory{side}=[linspace(trajectory{side}(1,1),trajectory{side}(2,1),50)',...
    linspace(trajectory{side}(1,2),trajectory{side}(2,2),50)',...
    linspace(trajectory{side}(1,3),trajectory{side}(2,3),50)'];

end



