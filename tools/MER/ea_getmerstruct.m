function [merstruct] = ea_getmerstruct(options)
%
% Returns merstruct from options
%
% Example: 
%   [merstruct] = ea_getmerstruct(options)

% __________________________________________________________________________________
% Copyright (C) 2017 University of Pittsburgh, Brain Modulation Lab
%
% Ari Kappel

if ~isfield(options,'native')
    options.native=0;
end
if ~isfield(options,'prefs')
    options.prefs=ea_prefs;
end
if ~isfield(options,'sides')
    options.sides=1:2;
end

uipatdirs=options.uipatdirs;

for pt=1:length(uipatdirs)
    options.uipatdirs{1}=uipatdirs{pt};
    [~,trajectory,markers]=ea_load_reconstruction(options);
    coords_mm=ea_resolvecoords(markers,options);
    merstruct.length = options.prefs.mer.length; %default is 24mm
    merstruct.offset = options.prefs.mer.offset; % default distance between mer tracts is 2mm
    merstruct.colormap = [0.5,0,0;0.5,0.5,0;0,0.5,0;0.5,0,0.5;0,0.5,0.5;0,0,0.5]; %Maroon,Olive,Green,Purple,Teal,Navy
    merstruct.defaultmer = ea_coords2mer(coords_mm,trajectory,merstruct,options);

for side = options.sides
    
    % contact_spacing = getfield(getappdata(resultfig,'elspec'),'contact_spacing');
    merstruct.defaultmer.central.trajectory{side} = ea_getmertrajectory(merstruct.defaultmer.central.coords_mm{side},0,merstruct.length,50);
    merstruct.defaultmer.anterior.trajectory{side} = ea_getmertrajectory(merstruct.defaultmer.anterior.coords_mm{side},0,merstruct.length,50);
    merstruct.defaultmer.posterior.trajectory{side} = ea_getmertrajectory(merstruct.defaultmer.posterior.coords_mm{side},0,merstruct.length,50);
    merstruct.defaultmer.lateral.trajectory{side} = ea_getmertrajectory(merstruct.defaultmer.lateral.coords_mm{side},0,merstruct.length,50);
    merstruct.defaultmer.medial.trajectory{side} = ea_getmertrajectory(merstruct.defaultmer.medial.coords_mm{side},0,merstruct.length,50);
    
end

end
function outputtrajectory = ea_getmertrajectory(coords_mm,dist,length,n)

if size(coords_mm,1)<2
    error('Check input a vector')
end
%%
dxyz = sqrt((diff(coords_mm(1:2,1))^2)+(diff(coords_mm(1:2,2))^2)+diff(coords_mm(1:2,3))^2);
slope = mean(diff(coords_mm))/dxyz;
startpoint = coords_mm(1,:)+slope.*dist;

outputtrajectory(:,1) = linspace(startpoint(1,1),startpoint(1,1)+slope(1)*length,n);
outputtrajectory(:,2) = linspace(startpoint(1,2),startpoint(1,2)+slope(2)*length,n);
outputtrajectory(:,3) = linspace(startpoint(1,3),startpoint(1,3)+slope(3)*length,n);


%% Equivalent solution
% __________________________________________________________________________________
%     normtrajvector=slope/norm(slope);
%     orth=null(normtrajvector);
%
%     startpoint.x = trajin(1,:)+orth(:,1)';
%     startpoint.y = trajin(1,:)+orth(:,2)';
%     startpoint=slope(1,:)-(2*(trajin(1,:)-slope(1,:)))
% __________________________________________________________________________________

function mer = ea_coords2mer(coords_mm,trajectory,merstruct,options)
        
offset = merstruct.offset;
contact_length = options.elspec.contact_length;

% x-axis --> negative = medial, positive = lateral
% y-axis --> negative = posterior, positive = anterior
% z-axis --> negative = inferior, positive = superior
for side=options.sides
    dxyz = pdist(coords_mm{side}(1:2,:),'euclidean');   %dxyz = sqrt((diff(coords_mm{side}(1:2,1))^2)+(diff(coords_mm{side}(1:2,2))^2)+diff(coords_mm{side}(1:2,3))^2);
    slope = mean(diff(coords_mm{side}))/dxyz;           %mean(diff(coords_mm{side}))/norm(mean(diff(coords_mm{side})))
    coords_mm{side} = coords_mm{side}-repmat(slope*contact_length/2,length(coords_mm{side}),1);
    trajectory{side} = trajectory{side}-repmat(slope*contact_length/2,length(trajectory{side}),1);
    if ~options.native
        mer.central.coords_mm{side} = coords_mm{side};
        mer.central.trajectory{side} = trajectory{side};
        mer.anterior.coords_mm{side} = [coords_mm{side}(:,1),coords_mm{side}(:,2)+offset,coords_mm{side}(:,3)]; %2mm anterior
        mer.anterior.trajectory{side} = [trajectory{side}(:,1),trajectory{side}(:,2)+offset,trajectory{side}(:,3)]; %2mm anterior
        mer.posterior.coords_mm{side} = [coords_mm{side}(:,1),coords_mm{side}(:,2)-offset,coords_mm{side}(:,3)]; %2mm posterior
        mer.posterior.trajectory{side} = [trajectory{side}(:,1),trajectory{side}(:,2)-offset,trajectory{side}(:,3)]; %2mm posterior
        if side==1 %right
            mer.lateral.coords_mm{side} = [coords_mm{side}(:,1)+offset,coords_mm{side}(:,2),coords_mm{side}(:,3)]; %2mm lateral (right is positive)
            mer.lateral.trajectory{side} = [trajectory{side}(:,1)+offset,trajectory{side}(:,2),trajectory{side}(:,3)]; %2mm lateral (right is positive)
            mer.medial.coords_mm{side} = [coords_mm{side}(:,1)-offset,coords_mm{side}(:,2),coords_mm{side}(:,3)]; %2mm medial (right is positive)
            mer.medial.trajectory{side} = [trajectory{side}(:,1)-offset,trajectory{side}(:,2),trajectory{side}(:,3)]; %2mm medial (right is positive)
        elseif side==2 %left
            mer.lateral.coords_mm{side} = [coords_mm{side}(:,1)-offset,coords_mm{side}(:,2),coords_mm{side}(:,3)]; %2mm lateral (left is negative)
            mer.lateral.trajectory{side} = [trajectory{side}(:,1)-offset,trajectory{side}(:,2),trajectory{side}(:,3)]; %2mm lateral (left is negative)
            mer.medial.coords_mm{side} = [coords_mm{side}(:,1)+offset,coords_mm{side}(:,2),coords_mm{side}(:,3)]; %2mm medial (left is negative)
            mer.medial.trajectory{side} = [trajectory{side}(:,1)+offset,trajectory{side}(:,2),trajectory{side}(:,3)]; %2mm medial (left is negative)
        end
    
    elseif options.native
        mer.central.coords_mm{side} = coords_mm{side};
        mer.central.trajectory{side} = trajectory{side};
        mer.anterior.coords_mm{side} = [coords_mm{side}(:,1),coords_mm{side}(:,2)+offset,coords_mm{side}(:,3)]; %2mm anterior
        mer.anterior.trajectory{side} = [trajectory{side}(:,1),trajectory{side}(:,2)+offset,trajectory{side}(:,3)]; %2mm anterior
        mer.posterior.coords_mm{side} = [coords_mm{side}(:,1),coords_mm{side}(:,2)-offset,coords_mm{side}(:,3)]; %2mm posterior
        mer.posterior.trajectory{side} = [trajectory{side}(:,1),trajectory{side}(:,2)-offset,trajectory{side}(:,3)]; %2mm posterior
        if side==1 %right
            mer.lateral.coords_mm{side} = [coords_mm{side}(:,1)+offset,coords_mm{side}(:,2),coords_mm{side}(:,3)]; %2mm lateral (right is positive)
            mer.lateral.trajectory{side} = [trajectory{side}(:,1)+offset,trajectory{side}(:,2),trajectory{side}(:,3)]; %2mm lateral (right is positive)
            mer.medial.coords_mm{side} = [coords_mm{side}(:,1)-offset,coords_mm{side}(:,2),coords_mm{side}(:,3)]; %2mm medial (right is positive)
            mer.medial.trajectory{side} = [trajectory{side}(:,1)-offset,trajectory{side}(:,2),trajectory{side}(:,3)]; %2mm medial (right is positive)
        elseif side==2 %left
            mer.lateral.coords_mm{side} = [coords_mm{side}(:,1)-offset,coords_mm{side}(:,2),coords_mm{side}(:,3)]; %2mm lateral (left is negative)
            mer.lateral.trajectory{side} = [trajectory{side}(:,1)-offset,trajectory{side}(:,2),trajectory{side}(:,3)]; %2mm lateral (left is negative)
            mer.medial.coords_mm{side} = [coords_mm{side}(:,1)+offset,coords_mm{side}(:,2),coords_mm{side}(:,3)]; %2mm medial (left is negative)
            mer.medial.trajectory{side} = [trajectory{side}(:,1)+offset,trajectory{side}(:,2),trajectory{side}(:,3)]; %2mm medial (left is negative)
        end
    end
end