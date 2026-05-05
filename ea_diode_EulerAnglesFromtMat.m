function angles = ea_diode_EulerAnglesFromtMat(T)
% Extract Euler angles from 4x4 affine
% Convention:
%   R = Rx(pitch) * Ry(yaw) * Rz(roll)

    A = T(1:3,1:3);

    % closest pure rotation
    [U,~,V] = svd(A);
    R = U*V';

    % fix reflections
    if det(R) < 0
        U(:,3) = -U(:,3);
        R = U*V';
    end

    % extract Euler angles
    angles.yaw   = asind(R(1,3));
    angles.pitch = atan2d(R(2,3), R(3,3));
    angles.roll  = atan2d(R(1,2), R(1,1));

    angles.roll = mod(angles.roll+180,360)-180;
    angles.R = R;
end