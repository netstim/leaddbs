function dir_angles = ea_diode_darkstar(roll,pitch,yaw,dirlevel,radius)
vizz = 0;
%% create vectors symbolizing the gaps between directional contacts at 60, 180 and 300 degrees
% and transform them to match lead trajectory and directional level
dirlevel = dirlevel(1:3);

ven = [0 0.65 -0.75 ]';
dor = [0 0.65 0.75 ]';
[M,~,~,~] = ea_diode_rollpitchyaw(roll-((2*pi)/6),pitch,yaw);
ven60 = M * ven;
dor60 = M * dor;
[M,~,~,~] = ea_diode_rollpitchyaw(roll-(3*(2*pi)/6),pitch,yaw);
ven180 = M * ven;
dor180 = M * dor;
[M,~,~,~] = ea_diode_rollpitchyaw(roll-(5*(2*pi)/6),pitch,yaw);
ven300 = M * ven;
dor300 = M * dor;

%% calculate intersecting points between vec60/180/300 and the z-plane through the dir-level artifact
vec60 = (dor60-ven60) / norm(dor60-ven60);              % unitvector from ven60 to dor60
dir_ven60 = dirlevel + ven60;                          % ventral point at 60° from the directional level
dir_dor60 = dirlevel + dor60;                          % dorsal point at 60° from the directional level
dir_x60 = (dirlevel(3) - dir_ven60(3)) / vec60(3);      % factor x of how many unitvectors dir_ven60 is distanced from the dirlevel in the z-dimension
dir_60 = dir_ven60 + (dir_x60 .* vec60);                % intersecting point of the line from ven60 to dor60 with the dirlevel plane in the z-dimension

vec180 = (dor180-ven180) / norm(dor180-ven180);
dir_ven180 = dirlevel + ven180;
dir_dor180 = dirlevel + dor180;
dir_x180 = (dirlevel(3) - dir_ven180(3)) / vec180(3);
dir_180 = dir_ven180 + (dir_x180 .* vec180);

vec300 = (dor300-ven300) / norm(dor300-ven300);
dir_ven300 = dirlevel + ven300;
dir_dor300 = dirlevel + dor300;
dir_x300 = (dirlevel(3) - dir_ven300(3)) / vec300(3);
dir_300 = dir_ven300 + (dir_x300 .* vec300);

% bla = [dir_60'; dir_180'; dir_300']
%% new stuff
      p1 = dir_60(1:2) - dirlevel(1:2);
      p2 = dir_180(1:2) - dirlevel(1:2);
      [dir_angles_new(1),dir_angles_new(2)] = calculatestreaks(p1,p2,radius);
      p1 = dir_180(1:2) - dirlevel(1:2);
      p2 = dir_300(1:2) - dirlevel(1:2);
      [dir_angles_new(3),dir_angles_new(4)] = calculatestreaks(p1,p2,radius);
      p1 = dir_300(1:2) - dirlevel(1:2);
      p2 = dir_60(1:2) - dirlevel(1:2);
      [dir_angles_new(5),dir_angles_new(6)] = calculatestreaks(p1,p2,radius);
      dir_angles_new = sort(dir_angles_new);
      
      dir_angles = dir_angles_new;
%%
% a=1+slope.^2;
% b=2*(slope.*(intercpt-centery)-centerx);
% c=centery.^2+centerx.^2+intercpt.^2-2*centery.*intercpt-radius.^2;
% 
% x=roots([a,b,c])';
% 
% %  Make NaN's if they don't intersect.
% 
% if ~isreal(x)
%     x=[NaN NaN]; y=[NaN NaN];
% else
%     y=[intercpt intercpt]+[slope slope].*x;
% end
% 
% % vertical slope case
% elseif abs(centerx-intercpt)>radius  % They don't intercept
%     x=[NaN;NaN]; y=[NaN;NaN];
%     else
%         x=[intercpt intercpt];
%         step=sqrt(radius^2-(intercpt-centerx)^2);
%         y=centery+[step,-step];
% end

%% create vectors corresponding to the dark lines of the artiface
dir_vec1 = (((dir_60 - dir_180) / norm(dir_60 - dir_180)) .* radius) + ((dir_60+dir_180) ./2) - dirlevel;
dir_vec1 = dir_vec1 ./ norm(dir_vec1);
dir_vec2 = (((dir_60 - dir_300) / norm(dir_60 - dir_300)) .* radius) + ((dir_60+dir_300) ./2) - dirlevel;
dir_vec2 = dir_vec2 ./ norm(dir_vec2);
dir_vec3 = (((dir_180 - dir_300) / norm(dir_180 - dir_300)) .* radius) + ((dir_180+dir_300) ./2) - dirlevel;
dir_vec3 = dir_vec3 ./ norm(dir_vec3);

% dir_vec1 = (dir_60 - dir_180) / norm(dir_60 - dir_180);
% dir_vec2 = (dir_60 - dir_300)  / norm(dir_60 - dir_300);
% dir_vec3 = (dir_180 - dir_300)  / norm(dir_180 - dir_300);

%% calculate the angles of the dark lines with respect to the y-axis
% dir_angles(1) = atan2(norm(cross(dir_vec1,[0 1 0])),dot(dir_vec1,[0 1 0]));
% if dir_vec1(1) < 0
%     dir_angles(1) = -dir_angles(1);
% end
% dir_angles(2) = atan2(norm(cross(dir_vec2,[0 1 0])),dot(dir_vec2,[0 1 0]));
% if dir_vec2(1) < 0
%     dir_angles(2) = -dir_angles(2);
% end
% dir_angles(3) = atan2(norm(cross(dir_vec3,[0 1 0])),dot(dir_vec3,[0 1 0]));
% if dir_vec3(1) < 0
%     dir_angles(3) = -dir_angles(3);
% end
% 
% dir_angles = [dir_angles (dir_angles + pi)];
% dir_angles(find(dir_angles>2*pi)) = dir_angles(find(dir_angles>2*pi)) - (2* pi);
% dir_angles(find(dir_angles<0)) = dir_angles(find(dir_angles<0)) + (2* pi);
% dir_angles = (2 *pi) -dir_angles;
% dir_angles = sort(dir_angles);

%% vizz

if vizz == 1
    figure
    temp = [dir_ven60'; dir_ven180'; dir_ven300'; dir_dor60'; dir_dor180'; dir_dor300'];
    scatter3(temp(:,1),temp(:,2),temp(:,3),'b')
    hold on
    plot3(temp([1 2 3 1],1),temp([1 2 3 1],2),temp([1 2 3 1],3),'b')
    plot3(temp([4 5 6 4],1),temp([4 5 6 4],2),temp([4 5 6 4],3),'b')
    clear temp
    
    axis equal
    temp = [dir_60'; dir_180'; dir_300'];
    scatter3(temp(:,1),temp(:,2),temp(:,3),'r')
    clear temp
    temp = [dir_60'; (dir_60 - dir_vec1)'; (dir_60 - dir_vec2)'; (dir_180 - dir_vec3)'];
    plot3(temp([1 2 3 1],1),temp([1 2 3 1],2),temp([1 2 3 1],3),'r')
    clear temp
    scatter3(dirlevel(1),dirlevel(2),dirlevel(3),'g');
    xlabel('x')
    ylabel('y')
    zlabel('z')
end
end

function [ws1,ws2] = calculatestreaks(p1,p2,radius)
      a = (p2(1) - p1(1))^2 + (p2(2) - p1(2))^2;
      b = 2 * (p1(1) * (p2(1) - p1(1)) + p1(2) * (p2(2) - p1(2)));
      c = p1(1) * p1(1) + p1(2) * p1(2) - radius^2;
      lambda1 = (-b + sqrt(b*b - 4*a*c)) / (2*a);
      lambda2 = (-b - sqrt(b*b - 4*a*c)) / (2*a);

      x1 = p1(1) + lambda1 * (p2(1) - p1(1));  %intersection of dark streak with validation circle
      y1 = p1(2) + lambda1 * (p2(2) - p1(2));
      x2 = p1(1) + lambda2 * (p2(1) - p1(1));
      y2 = p1(2) + lambda2 * (p2(2) - p1(2));

%       ws1 = rad2deg(atan2(y1, x1));
%       ws2 = rad2deg(atan2(y2, x2));
%       ws1 = -90 + ws1;            % angle clockwise with respect to +y
%       ws2 = -90 + ws2;
%       if ws1 < 0
%           ws1 = ws1 + 360;
%       end
%       if ws2 < 0 
%           ws2 = ws2 + 360;
%       end
     ws1 = atan2(y1, x1);
      ws2 = atan2(y2, x2);
      ws1 = -(pi/2) + ws1;            % angle clockwise with respect to +y
      ws2 = -(pi/2) + ws2;
      if ws1 < 0
          ws1 = ws1 + (2*pi);
      end
      if ws2 < 0 
          ws2 = ws2 + (2*pi);
      end
end