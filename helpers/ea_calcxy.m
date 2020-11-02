function [xunitv, yunitv] = ea_calcxy(head, tail, y)
% Calculate the unit vectors of x and y markers of the lead
%
% head and tail coordinates should be accurate. The rotation/orientation of
% the y marker in the x-y plane is more important than the absolute
% coordinate. It will be recalculated/adjusted.
%
% If y is not available, can choose whether to force y axis pointing
% anterior, or to calculate x and y axes using null function (orthonormal
% basis of the null space of the z axis).

if ~exist('y', 'var')
    y = 'anterior';
end

trajvector = diff([head; tail]); % z axis, poiting to the top of the lead
normtrajvector = trajvector/norm(trajvector); % Unit vector along z axis
if isnumeric(y)
    x = cross(y/norm(y),[0 0 1]); % [1 0 0] rotated according to y
    x = x - dot(x,normtrajvector) * normtrajvector; % Project x down to the trajectory
    xunitv = x/norm(x); % Calculate unit vector along x axis
    yunitv = -cross(x,normtrajvector); % Calculate unit vector along y axis
elseif ischar(y)
    switch y
        case 'anterior' % Force y axis pointing anterior
            y = [0, normtrajvector(3), -normtrajvector(2)];
            x = cross(y, normtrajvector);
            yunitv = y/norm(y);
            xunitv = x/norm(x);
        case 'null' % Use orthonormal basis calculated from trajvector
            ortho = null(normtrajvector)';
            yunitv = ortho(1,:)/norm(ortho(1,:));
            xunitv = -ortho(2,:)/norm(ortho(2,:));
    end
end
