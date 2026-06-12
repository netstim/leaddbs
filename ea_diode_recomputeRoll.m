function angles = ea_diode_recomputeRoll(head,tail,y)
    % by Till Dembek & ChatGPT 5/2026
    yvec = y - head;
    R = leadFrameFromVectors(head, tail, yvec);
    % orthogonalize
    [U,~,V] = svd(R);
    R = U*V';

      % trajectory angles
    % angles.yaw   = asind(R(1,3));           % +right entry
    % angles.pitch = atan2d(R(2,3), R(3,3));  % +anterior entry    
    angles.deg.roll  = atan2d(-R(1,2), R(1,1));  % roll: CW positive viewed from inferior

    % normalize
    angles.deg.roll = mod(angles.deg.roll+180,360)-180;
    angles.rad.roll = deg2rad(angles.deg.roll);
end

function R = leadFrameFromVectors(head, tail, yvec)
    % z-axis = lead axis (superior direction)
    z = (tail(:) - head(:));
    z = z / norm(z);

    % y-axis = orientation direction
    y = yvec(:);
    y = y / norm(y);

    % ensure orthogonal projection
    y = y - dot(y,z)*z;
    y = y / norm(y);

    % x-axis via right-handed frame
    x = cross(y, z);
    x = x / norm(x);

    % recompute y for exact orthogonality
    y = cross(z, x);

    % rotation matrix
    R = [x y z];
end