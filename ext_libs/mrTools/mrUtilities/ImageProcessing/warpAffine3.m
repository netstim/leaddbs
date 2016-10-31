function result = warpAffine3(in,A,badVal,B,method,outSize)
%
% function result = warpAffine3(in,A,[badVal],[border],[method],[outSize])
%
% in: input volume, 3D array
% A: 3x4 affine transform matrix or a 4x4 matrix with [0 0 0 1]
%    for the last row.
% badVal: if a transformed point is outside of the volume, badVal is used
%    (default = 0).
% border: number of voxels to put in the border (default = 0, no border).
% method:  'nearest', 'linear', 'cubic', or 'spline' (default = 'linear')
% outSize: optional specification of output size defaults to input
%          size. Should be a 1x3 vector like what size returns
%
% result: output volume
%
%

if ieNotDefined('badVal')
    badVal = 0;
end
if ieNotDefined('B')
    B = 0;
end
if ieNotDefined('method')
    method = 'linear';
end
if ieNotDefined('outSize')
  outSize = size(in);
end
if (size(A,1)>3)
    A = A(1:3,:);
end

% original size
NyO = outSize(1);
NxO = outSize(2);
NzO = outSize(3);

% Compute coordinates corresponding to input volume
% and transformed coordinates for result
[xgrid,ygrid,zgrid]=meshgrid(1:NxO,1:NyO,1:NzO);
coords=[xgrid(:)'; ygrid(:)'; zgrid(:)'];
homogeneousCoords=[coords; ones(1,size(coords,2))];
warpedCoords=A*homogeneousCoords;

% Interpolate
if B~=0
    % Put a border by replicating edge voxels.
    for k =1:size(in,3)
        inB(:,:,k) = addBorder(in(:,:,k),B,B,3);
    end
    % Repeat the first and last slices
    in = cat(3,repmat(inB(:,:,1),[1 1 B]), inB, repmat(inB(:,:,end),[1 1 B]));
    % Coordinates corresponding to the volume with borders.
    [xB,yB,zB]=meshgrid(1-B:size(in,2)-B,1-B:size(in,1)-B,1-B:size(in,3)-B);
    % Use interp3 for interpolation
    result = interp3(xB,yB,zB,in,warpedCoords(1,:),warpedCoords(2,:),warpedCoords(3,:),...
        method, badVal);
else
    if strcmp(method,'linear')
        % If no border added and linear interpolation then use myCinterp3
        if isreal(in)
            result = myCinterp3(double(in), [size(in,1) size(in,2)], size(in,3), ...
                [warpedCoords(1,:);warpedCoords(2,:);warpedCoords(3,:)]', badVal);
        else
            inReal = real(in);
            inImag = imag(in);
            resultReal = myCinterp3(double(inReal), [size(in,1) size(in,2)], size(in,3), ...
                [warpedCoords(1,:);warpedCoords(2,:);warpedCoords(3,:)]', badVal);
            resultImag = myCinterp3(double(inImag), [size(in,1) size(in,2)], size(in,3), ...
                [warpedCoords(1,:);warpedCoords(2,:);warpedCoords(3,:)]', badVal);
            result = resultReal + j*resultImag;
        end
    else
        result = interp3(in,warpedCoords(1,:),warpedCoords(2,:),warpedCoords(3,:),...
            method, badVal);
    end
end

% Reshape result
result = reshape(result,[NyO NxO NzO]);

return;

%%% Debug

slice=[1 2 3; 4 5 6; 7 8 9]';
slice=[1 1 1; 3 3 3; 5 5 5]';
input=ones(3,3,4);
for z=1:4
    input(:,:,z)=slice;
end

A= [1 0 0 .5;
    0 1 0 0;
    0 0 1 0;
    0 0 0 1];

A= [1 0 0 .5;
    0 1 0 .5;
    0 0 1 0];

res=warpAffine3(input,A)
res=warpAffine3(input,A,-100)
res=warpAffine3(input,A,NaN)
resB=warpAffine3(input,A,-100,1)

for z=1:4
    input(:,:,z)=z*ones(3,3);
end
A= [1 0 0 0;
    0 1 0 0;
    0 0 1 .5;
    0 0 0 1];
res=warpAffine3(input,A)
resB=warpAffine3B(input,A,NaN,1)

theta=.02;
A=[cos(theta) sin(theta) 0 0;
    -sin(theta) cos(theta) 0 0;
    0 0 1 0;
    0 0 0 1];
input=rand(5,5,5);
res1a=warpAffine3(input,A)
res2a=warpAffine3(input,A,-100)
res3a=warpAffine3(input,A,-100,1)

