function angles = ea_diode_trajectoryToYawPitchPolar(head,tail)
% by Till Dembek & ChatGPT 5/2026
% Compute trajectory angles for DBS leads
%
% Assumptions:
%   - RAS coordinate system
%   - head = inferior tip
%   - tail = superior point on lead
%   - trajectory vector = tail - head
%
% Sign conventions:
%   +pitch = anterior entry point
%   +yaw   = rightward entry point

% normalize
v = tail-head;
v = v(:) / norm(v);

vx = v(1);   % right
vy = v(2);   % anterior
vz = v(3);   % superior

% sanity check
if vz < 0
    warning('Trajectory points inferiorly; expected tail-head to point superiorly.');
end

% angle calculations deg
angles.deg.yaw   = asind(vx);        % + rightward entry
angles.deg.pitch = atan2d(vy, vz);   % + anterior entry
angles.deg.polar = atan2d(sqrt(vx^2 + vy^2), vz); % total tilt from z-axis
% angle calculations rad
angles.rad.yaw   = asin(vx);        % + rightward entry
angles.rad.pitch = atan2(vy, vz);   % + anterior entry
angles.rad.polar = atan2(sqrt(vx^2 + vy^2), vz); % total tilt from z-axis
% also return unitvector
angles.unitvector = v;
end