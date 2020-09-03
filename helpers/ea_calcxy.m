function [x, y] = ea_calcxy(head, tail, y)
% Calculate the unit vectors of x and y markers of the lead
%
% head and tail coordinates should be accurate. The rotation/orientation of
% the y marker in the x-y plane is more important than the absolute
% coordinate. It will be recalculated/adjusted.

trajvector = diff([head; tail]); % z axis, poiting to the top of the lead
normtrajvector = trajvector/norm(trajvector); % Unit vector along z axis
y = y - dot(y,normtrajvector)*normtrajvector; % Adjust y axis
x = cross(y, normtrajvector); % Calculate x axis based on y and z axis
y = y/norm(y); % Calculate unit vector along y axis
x = x/norm(x); % Calculate unit vector along x axis
