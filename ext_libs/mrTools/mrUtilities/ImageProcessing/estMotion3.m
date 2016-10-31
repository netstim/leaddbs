function [M,w3d] = estMotion3(vol1,vol2,rotFlag,robustFlag,phaseFlag,crop,CB,SC,mask)
%
% function [M,w] = estMotion3(vol1,vol2,[rotFlag],[robustFlag],[phaseFlag],[crop],[CB],[SC])
%
% vol1 and vol2 are volumes, 3d arrays which can be either real or
% complex-valued
%
% Solves fs^t theta + ft = 0
% where theta = B p is image velocity at each pixel
%       B is 3x6 (3x12 if affine) matrix that depends on image positions
%       p is vector of trans+rot (or affine) motion parameters
%       fs is vector of spatial derivatives at each pixel
%       ft is temporal derivative at each pixel
% Multiplying fs^t B gives a 1x6 (1x12 if affine) vector for each pixel.  Piling
% these on top of one another gives A, an Nx6 (Nx12 if affine) matrix, where N is
% the number of pixels. Solves M p = -ftVol where ftVol is an Nx1
% vector of the the temporal derivatives at every pixel.
%
% M is 4x4 (rotation+translation) transform matrix: X' = M X
% where X=(x,y,z,1) is starting position in homogeneous coords
% and X'=(x',y',z',1) is ending position
%
% rotFlag: If True, then M is a rotation+translation, otherwise, is a
% general affine transform. Default: 1.
%
% robustFlag: If True, uses robust M-estimator (see robustMest) with
% parameters CB and SC. Default: 0.
%
% w is a volume (same size as vol1 and vol2) containing the weights
% returned by the robust estimator (see robustMest.m).
%
% phaseFlag: If True, treats the input volumes as containing the phase of
% complex-values (on the unit circle in the complex plane) and computes the
% motion estimates on the complex values. This is useful if the image
% intensities in the volumes are from a periodic domain such that the
% "brightest" and "darkest" intensities are actually very similar to each
% other. Default: 0.
%
% crop specifies border size to crop/ignore  around all sides of the volume.
%     Should be of the form [ymin xmin zmin; ymax xmax zmax]
%     Default crops 2 pixel border: [2 2 2; (size(vol1) - [2 2 2])].
%
% ON 08/00 - modified to permit registration of complex valued volumes

% default values
if ~exist('rotFlag','var')
    rotFlag = 1;
end
if ~exist('robustFlag','var')
    robustFlag = 0;
end
if ~exist('phaseFlag','var')
    phaseFlag = 0;
end
if ~exist('crop','var') || isempty(crop)
    crop = [1 1 1; size(vol1)];
end
if ~exist('CB','var')
    CB = [];
end
if ~exist('SC','var')
    SC = [];
end
if ~exist('mask','var') 
    mask = [];
end

% Phase
if phaseFlag
    vol1 = exp(j*vol1);
    vol2 = exp(j*vol2);
end

% Compute derivatives
%[fxVol,fyVol,fzVol,ftVol] = computeDerivatives3(vol1,vol2,'same'); %this is a departure from what was done before ('valid' was used), just trying it out
%And it is necessary in order to use the mask without having to figure out the dimensions to pad with NaNs
%although we know that it's 2 on each side in all dimensions because computeDerivatives3 uses kernels of length 5, so:
[fxVol,fyVol,fzVol,ftVol] = computeDerivatives3(vol1,vol2,'valid'); 

% *** We assume 5tap derivative filters.
filterTaps = 5;
nInvalidVoxels = floor(filterTaps/2);
if ~isempty(mask) 
   %reduce the mask to account fo the use of the 'valid' derivation
   mask = mask(3:end-nInvalidVoxels,3:end-nInvalidVoxels,3:end-nInvalidVoxels);
   %apply the mask
   fxVol(~mask)=NaN;
   fyVol(~mask)=NaN;
   fzVol(~mask)=NaN;
   ftVol(~mask)=NaN;
end

% Adjust crop
crop(1,:) = max(crop(1,:)-nInvalidVoxels,[1 1 1]);
crop(2,:) = min(crop(2,:)-nInvalidVoxels,size(vol1)-2*nInvalidVoxels);

% Meshgrid on original volume. Then below we throw the edges of the
% meshgrid to make it the same size as the derivative volumes.
[xgrid,ygrid,zgrid] = meshgrid(1:size(vol1,2),1:size(vol1,1),1:size(vol1,3));
xgrid = xgrid(3:end-nInvalidVoxels,3:end-nInvalidVoxels,3:end-nInvalidVoxels);
ygrid = ygrid(3:end-nInvalidVoxels,3:end-nInvalidVoxels,3:end-nInvalidVoxels);
zgrid = zgrid(3:end-nInvalidVoxels,3:end-nInvalidVoxels,3:end-nInvalidVoxels);


% Subsample and crop
indicesY = crop(1,1):2:crop(2,1);
indicesX = crop(1,2):2:crop(2,2);
indicesZ = crop(1,3):2:crop(2,3);
fxVol = fxVol(indicesY,indicesX,indicesZ);
fyVol = fyVol(indicesY,indicesX,indicesZ);
fzVol = fzVol(indicesY,indicesX,indicesZ);
ftVol = ftVol(indicesY,indicesX,indicesZ);
xgrid = xgrid(indicesY,indicesX,indicesZ);
ygrid = ygrid(indicesY,indicesX,indicesZ);
zgrid = zgrid(indicesY,indicesX,indicesZ);

dimsS=size(fxVol);
pts=find((~isnan(fxVol))&(~isnan(fyVol))&(~isnan(fzVol))&(~isnan(ftVol)));
%disp(['numPts=',num2str(length(pts))]);
fxVol = fxVol(pts);
fyVol = fyVol(pts);
fzVol = fzVol(pts);
ftVol = ftVol(pts);
xVol = xgrid(pts);
yVol = ygrid(pts);
zVol = zgrid(pts);

if rotFlag
    A = [fxVol(:),  fyVol(:), fzVol(:),...
        fzVol(:).*yVol(:)-fyVol(:).*zVol(:) ...
        fxVol(:).*zVol(:)-fzVol(:).*xVol(:)...
        fyVol(:).*xVol(:)-fxVol(:).*yVol(:)];
else
    A = [xVol(:).*fxVol(:), yVol(:).*fxVol(:), zVol(:).*fxVol(:), fxVol(:),...
        xVol(:).*fyVol(:), yVol(:).*fyVol(:), zVol(:).*fyVol(:), fyVol(:),...
        xVol(:).*fzVol(:), yVol(:).*fzVol(:), zVol(:).*fzVol(:), fzVol(:)];
end

b = -ftVol(:);

% This allows to solve for complex valued volumes because the motion
% parameters must be real. It is equivalent to consider that the imaginary
% part is adding new constraints.
complexFlag = (~isreal(A))|(~isreal(b));
if complexFlag
    A = [real(A); imag(A)];
    b = [real(b); imag(b)];
end

if robustFlag
    [p w] = robustMest(A,b,CB,SC);
    if complexFlag
        % rearrange the weights and return them in complex form
        w = w(1:length(w)/2) + sqrt(-1) * w(length(w)/2+1:end);
    end
    w3d = zeros(dimsS);
    w3d(pts) = w;
else
    p = A\b;
    w3d = [];
end

if rotFlag
    M = [quatR2mat(quatrot(p(4:6))) p(1:3); 0 0 0 1];
else
    M= [1+p(1)  p(2)   p(3)   p(4);
        p(5)   1+p(6)  p(7)   p(8);
        p(9)    p(10) 1+p(11) p(12);
        0       0     0      1];
end

return;

%%%%%%%%%
% Debug %
%%%%%%%%%

% test with translation
in=rand(30,32,28);
A= [1 0 0 .2;
    0 1 0 .3;
    0 0 1 .4;
    0 0 0 1];
vol1=warpAffine3(in,A);
vol2=warpAffine3(in,inv(A));
A*A
crop = [2 2 2; (size(vol1) - [2 2 2])];
% default - rot and LS
estMotion3(vol1,vol2,0,0,0,crop)
% rot and robust
estMotion3(vol1,vol2,1,1,0,crop)
% affine and LS
estMotion3(vol1,vol2,0,0,0,crop)
% affine and robust
estMotion3(vol1,vol2,0,1,0,crop)

% test with rotation
theta=.02;
A=[cos(theta) sin(theta) 0 0;
    -sin(theta) cos(theta) 0 0;
    0 0 1 0;
    0 0 0 1];
vol1=warpAffine3(in,A);
vol2=warpAffine3(in,inv(A));
A*A
% default - rot and LS
estMotion3(vol1,vol2,0,0,0,crop)
% rot and robust
estMotion3(vol1,vol2,1,1,0,crop)
% affine and LS
estMotion3(vol1,vol2,0,0,0,crop)
% affine and robust
estMotion3(vol1,vol2,0,1,0,crop)

% test complex valued volumes
in = rand(30,32,28) + j*rand(30,32,28);
theta = .02;
A = [cos(theta) sin(theta) 0 .2;
    -sin(theta) cos(theta) 0 .3;
     0 0 1 .4;
     0 0 0 1];
vol1 = warpAffine3(in,A);
vol2 = warpAffine3(in,inv(A));
A*A
crop = [4 4 4; (size(vol1) - [4 4 4])];
estMotion3(vol1,vol2,1,0,crop)
estMotionIter3(vol1,vol2,3,eye(4),1,0,0,crop)

% test phase-only volumes
phase = 2*pi*rand(30,32,28) - pi;
vol1 = warpAffine3(phase,A);
vol2 = warpAffine3(phase,inv(A)) + 2*pi;
A*A
crop = [4 4 4; (size(vol1) - [4 4 4])];
% these 2 should fail
estMotion3(vol1,vol2,1,0,0,crop)
estMotionIter3(vol1,vol2,3,eye(4),1,0,0,crop)
% these 2 should work
estMotion3(vol1,vol2,1,0,1,crop)
estMotionIter3(vol1,vol2,3,eye(4),1,0,1,crop)

%%%%%%%%%%%%%%%%%%%%%%
% test with outliers %
%%%%%%%%%%%%%%%%%%%%%%

% translation
dims=[30,32,40];
in=rand(dims);
A= [1 0 0 .2;
    0 1 0 .3;
    0 0 1 .4;
    0 0 0 1];
vol1 = warpAffine3(in,A);
vol2 = warpAffine3(in,inv(A));
% putting inconsistent information in upper left corner of vol2
Nc=3;
vol2(1:round(dims(1)/Nc), 1:round(dims(2)/Nc), :) = ...
    rand(round(dims(1)/Nc), round(dims(2)/Nc), dims(3));
% default - rot and LS
ArotLS = estMotion3(vol1,vol2);
% rot and robust
[ArotRob wr] = estMotion3(vol1,vol2,1,1);
% affine and LS
AaffLS = estMotion3(vol1,vol2,0,0);
% affine and robust
[AaffRob wa] = estMotion3(vol1,vol2,0,1);
ArotLS
ArotRob
AaffLS
AaffRob
dimw = size(wr);
imagesc([reshape(wr,[dimw(1) dimw(2)*dimw(3)]); reshape(wa,[dimw(1) dimw(2)*dimw(3)])]);
colormap(gray);axis('image');axis('off')

% test with rotation
dims = [30, 32, 28];
in=rand(dims);
theta=.02;
A=[cos(theta) sin(theta) 0 0;
    -sin(theta) cos(theta) 0 0;
    0 0 1 0;
    0 0 0 1];
vol1=warpAffine3(in,A);
vol2=warpAffine3(in,inv(A));
% putting inconsistent information in upper left corner of vol2
Nc=3;
vol2(1:round(dims(1)/Nc), 1:round(dims(2)/Nc), :) = ...
    rand(round(dims(1)/Nc), round(dims(2)/Nc), dims(3));
% default - rot and LS
ArotLS = estMotion3(vol1,vol2);
% rot and robust
[ArotRob wr] = estMotion3(vol1,vol2,1,1);
% affine and LS
AaffLS = estMotion3(vol1,vol2,0,0);
% affine and robust
[AaffRob wa] = estMotion3(vol1,vol2,0,1);
A*A
ArotLS
ArotRob
AaffLS
AaffRob
dimw = size(wr);
imagesc([reshape(wr,[dimw(1) dimw(2)*dimw(3)]); reshape(wa,[dimw(1) dimw(2)*dimw(3)])]);
colormap(gray);axis('image');axis('off')
