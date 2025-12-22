function R = LeG_genRotMat(vec1,vec2)

% Generates rotation matrix R for rotating vec1 to match vec2 in 3D space.

vec1 = vec1/sqrt(sum(vec1.^2)); %unit vector in start direction
vec2 = vec2/sqrt(sum(vec2.^2)); %unit vector in end direction
Theta = -acos(dot(vec1,vec2)); %angle between vec1 and vec2
RotAxis = cross(vec1,vec2); RotAxis = (RotAxis/sqrt(sum(RotAxis.^2))); %normalized axis of rotation
if any(isnan(RotAxis))
    R = eye(4);
else
    R = [cos(Theta)+RotAxis(1)^2*(1-cos(Theta)),RotAxis(1)*RotAxis(2)*(1-cos(Theta))-RotAxis(3)*sin(Theta),RotAxis(1)*RotAxis(3)*(1-cos(Theta))+RotAxis(2)*sin(Theta);
        RotAxis(2)*RotAxis(1)*(1-cos(Theta))+RotAxis(3)*sin(Theta),cos(Theta)+RotAxis(2)^2*(1-cos(Theta)),RotAxis(2)*RotAxis(3)*(1-cos(Theta))-RotAxis(1)*sin(Theta);
        RotAxis(3)*RotAxis(1)*(1-cos(Theta))-RotAxis(2)*sin(Theta),RotAxis(3)*RotAxis(2)*(1-cos(Theta))+RotAxis(1)*sin(Theta),cos(Theta)+RotAxis(3)^2*(1-cos(Theta))]; %Rotate about RotAxis by Theta
    
    R = [[R',zeros(3,1)];[zeros(1,3),1]];
end
