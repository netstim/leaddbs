function qr = quatrot(w);
% QUATROT - returns the quaternion representing the rotation given by
%           the 3D vector w (rotation of ||w|| around the direction given by w)
%
%   qr = quatrot(w);
%
%  ON - 3/99

th = norm(w);
% avoid divide by zero error
if (th == 0)
    qr = [1 0*w(:)'];
else  
    qr = [cos(th/2) sin(th/2)*w(:)'/th];
end;
