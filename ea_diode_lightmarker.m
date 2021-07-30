function dir_angles = ea_diode_lightmarker(roll,pitch,yaw,marker)
vizz = 0;
%% create vectors symbolizing the gaps between directional contacts at 60, 270 and 300 degrees
% and transform them to match lead trajectory and directional level
marker = marker(1:3);

ven = [0 0.65 -0.75 ]';
dor = [0 0.65 0.75 ]';

[M,~,~,~] = ea_diode_rollpitchyaw(roll-(pi/2),pitch,yaw);
ven90 = M * ven;
dor90 = M * dor;
[M,~,~,~] = ea_diode_rollpitchyaw(roll-(3*(pi/2)),pitch,yaw);
ven270 = M * ven;
dor270 = M * dor;
%% calculate intersecting points between vec60/270/300 and the z-plane through the dir-level artifact
vec90 = (dor90-ven90) / norm(dor90-ven90);              % unitvector from ven0 to dor0
dir_ven90 = marker + ven90;                          % ventral point at 0° from the directional level
dir_dor90 = marker + dor90;                          % dorsal point at 0° from the directional level
dir_x90 = (marker(3) - dir_ven90(3)) / vec90(3);      % factor x of how many unitvectors dir_ven0 is distanced from the marker in the z-dimension
dir_90 = dir_ven90 + (dir_x90 .* vec90);                % intersecting point of the line from ven0 to dor0 withe the marker plane in the z-dimension

vec270 = (dor270-ven270) / norm(dor270-ven270);
dir_ven270 = marker + ven270;
dir_dor270 = marker + dor270;
dir_x270 = (marker(3) - dir_ven270(3)) / vec270(3);
dir_270 = dir_ven270 + (dir_x270 .* vec270);

%% create vectors corresponding to the dark lines of the artiface
dir_vec1 = (dir_90 - dir_270) / norm(dir_90 - dir_270);

%% calculate the angles of the dark lines with respect to the y-axis
dir_angles(1) = atan2(norm(cross(dir_vec1,[0 1 0])),dot(dir_vec1,[0 1 0]));
if dir_vec1(1) < 0
    dir_angles(1) = -dir_angles(1);
end    

dir_angles = [dir_angles (dir_angles + pi)];
dir_angles(find(dir_angles>2*pi)) = dir_angles(find(dir_angles>2*pi)) - (2* pi);
dir_angles(find(dir_angles<0)) = dir_angles(find(dir_angles<0)) + (2* pi);
dir_angles = (2 *pi) -dir_angles;
dir_angles = sort(dir_angles);

%% vizz

if vizz == 1
    figure
    temp = [dir_ven60'; dir_ven270'; dir_ven300'; dir_dor60'; dir_dor270'; dir_dor300'];
    scatter3(temp(:,1),temp(:,2),temp(:,3),'b')
    hold on
    plot3(temp([1 2 3 1],1),temp([1 2 3 1],2),temp([1 2 3 1],3),'b')
    plot3(temp([4 5 6 4],1),temp([4 5 6 4],2),temp([4 5 6 4],3),'b')
    clear temp
    
    axis equal
    temp = [dir_60'; dir_270'; dir_300'];
    scatter3(temp(:,1),temp(:,2),temp(:,3),'r')
    clear temp
    temp = [dir_60'; (dir_60 - dir_vec1)'; (dir_60 - dir_vec2)'; (dir_270 - dir_vec3)'];
    plot3(temp([1 2 3 1],1),temp([1 2 3 1],2),temp([1 2 3 1],3),'r')
    clear temp
    scatter3(marker(1),marker(2),marker(3),'g');
    xlabel('x')
    ylabel('y')
    zlabel('z')
end
end