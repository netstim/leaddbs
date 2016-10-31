function xform = shiftOriginXform(origin)
%
% function xform = shiftOriginXform([origin])
%
% Returns 4x4 homogeneous transform matrix that shifts the coordinate frame
% to a different origin. This is used by a number of functions in mrLoadRet
% and mrAlign because matlab indexes from 1 but nifti uses 0,0,0 as the
% origin.
%
% origin: option 3 vector to specify the shift in origin (default [0,0,0]);
%
% djh 1/2007

if ieNotDefined('origin')
    origin = [0,0,0]';
end

xform = eye(4);
xform(1:3,4) = origin(:) - 1;
