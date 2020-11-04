function [rotation, yunitv] = ea_calc_rotation(y, head)
% Calc the rotation of y marker in the x-y plane

yvec = y - head;
yvec(3) = 0;
yunitv = yvec ./ norm(yvec);
rotation = rad2deg(atan2(norm(cross([0 1 0],yunitv)),dot([0 1 0],yunitv)));

% Rotation to right side is negative; Rotation to left side is positive
if y(1) > head(1)
    rotation = -rotation;
end
