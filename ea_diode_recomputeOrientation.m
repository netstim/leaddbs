function angles = ea_diode_recomputeOrientation(head, tail, y)
% Till Dembek & ChatGPT 05/2026
% Orientation angle relative to anterior (+y)
% 0° = anterior, +90° = right, CW positive from inferior
    yvec = y(:)-head(:);
    % lead axis
    z = tail(:) - head(:);
    z = z / norm(z);

    % projected anterior
    a = [0;1;0];
    a = a - dot(a,z)*z;
    a = a / norm(a);

    % projected orientation vector
    yvec = yvec - dot(yvec,z)*z;
    yvec = yvec / norm(yvec);

    % signed angle
    angles.deg.orientation = atan2d(dot(cross(a,yvec),z), dot(a,yvec));

    % optional: compass format
    angles.deg.orientation = mod(angles.deg.orientation+180,360)-180;
    angles.rad.orientation = deg2rad(angles.deg.orientation);
end