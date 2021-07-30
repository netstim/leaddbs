function [M,Mz,My,Mx] = ea_diode_rollpitchyaw(roll,pitch,yaw)
%% returns a 3x3 rotation matrix for roll->pitch->yaw rotation in that order
%% all in radians
a = pitch; %around x axis
b = yaw; %around y axis
c = roll; %around z axis

Mx = [1, 0, 0; 0, cos(a), sin(a); 0, -sin(a), cos(a)];
My = [cos(b), 0, sin(b); 0, 1, 0; -sin(b), 0, cos(b)];
Mz = [cos(c), -sin(c), 0; sin(c), cos(c), 0; 0, 0, 1];

M = Mx * My * Mz;
end