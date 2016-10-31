function volInterp = interpVolume(vol, xform, inplaneSize, badval)
% function volInterp = interpVolume(volume, xform, inplaneSize, badval)
%
%        $Id$
%
% djh, 1/2007
% jb 10/2010: clear unused variables to save memory usage

if ieNotDefined('badval')
    badval = NaN;
end

% Shift xform: matlab indexes from 1 but nifti uses 0,0,0 as the origin. 
shiftXform = shiftOriginXform;
xform = shiftXform \ xform * shiftXform;

% Generate coordinates
[y,x,z] = meshgrid(1:inplaneSize(2), 1:inplaneSize(1), 1:inplaneSize(3));
dims = inplaneSize;
numPixels = prod(dims);
xvec = reshape(x,1,numPixels);
yvec = reshape(y,1,numPixels);
zvec = reshape(z,1,numPixels);
clear('x','y','z');
coords = [xvec; yvec; zvec; ones(1,numPixels)];
clear('xvec','yvec','zvec');

% Transform coordinates
coordsXform = xform * coords;
clear('coords');
xi = reshape(coordsXform(1,:),dims);
yi = reshape(coordsXform(2,:),dims);
zi = reshape(coordsXform(3,:),dims);
clear('coordsXform');
% Interpolate
% Note: interp3 treats x and y in right-handed coordinate system, not in
% matrix index order so we need to swap them here. See example code below.
volInterp = interp3(vol,yi,xi,zi,'linear',badval);

return

% Example code for how to use interp3
bar = [1:27]
foo = reshape(bar,[3,3,3])
x = 2 * ones(3,3)
[z,y] = meshgrid(1:3,1:3)
squeeze(foo(2,:,:))
interp3(foo,y,x,z)

% Test/Debug
vol = foo;
vol = ALIGN.inplanes;
volInterp = interpVolume(vol,eye(4),size(vol),0);
difference = (vol - volInterp) ./ vol;
max(difference(:))
min(difference(:))
