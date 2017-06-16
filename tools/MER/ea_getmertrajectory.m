function outputtrajectory = ea_getmertrajectory(coords_mm,dist,length,n)

% __________________________________________________________________________________
% Copyright (C) 2017 University of Pittsburgh, Brain Modulation Lab
%
% Ari Kappel


if size(coords_mm,1)<2
    error('Must input a vector')
end
%%
dxyz = sqrt((diff(coords_mm(1:2,1))^2)+(diff(coords_mm(1:2,2))^2)+diff(coords_mm(1:2,3))^2);
slope = mean(diff(coords_mm))/dxyz;
startpoint = coords_mm(1,:)+slope.*dist;

%% Equivalent solution
%     normtrajvector=slope/norm(slope);
%     orth=null(normtrajvector);
%
%     startpoint.x = trajin(1,:)+orth(:,1)';
%     startpoint.y = trajin(1,:)+orth(:,2)';
%     startpoint=slope(1,:)-(2*(trajin(1,:)-slope(1,:)))

outputtrajectory(:,1) = linspace(startpoint(1,1),startpoint(1,1)+slope(1)*length,n);
outputtrajectory(:,2) = linspace(startpoint(1,2),startpoint(1,2)+slope(2)*length,n);
outputtrajectory(:,3) = linspace(startpoint(1,3),startpoint(1,3)+slope(3)*length,n);