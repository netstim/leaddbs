function [M,w2d] = estMotion2(im1,im2,rotFlag,robustFlag,crop,CB,SC)
%
% function [M,w] = estMotion2(im1,im2,rotFlag,robustFlag,[crop],[CB],[SC])
%
% im1 and im2 are images
%
% M is 3x3 transform matrix: X' = M X
% where X=(x,y,1) is starting position in homogeneous coords
% and X'=(x',y',1) is ending position
%
% Solves fs^t theta + ft = 0
% where theta = B p is image velocity at each pixel
%       B is 2x3 (2x6 if affine) matrix that depends on image positions
%       p is vector of trans+rot (or affine) motion parameters
%       fs is vector of spatial derivatives at each pixel
%       ft is temporal derivative at each pixel
% Mulitplying fs^t B gives a 1x3 (1x6 if affine) vector for each pixel.  Piling
% these on top of one another gives A, an Nx3 (Nx6 if affine) matrix, where N is
% the number of pixels.  Solve M p = ft where ft is now an Nx1
% vector of the the temporal derivatives at every pixel.
%
% If rotFlag is activated (~=0), then M is a rotation+translation,
% otherwise, is a general affine transform
%
% If robustFlag is activated (~=0) then uses a robust M-estimator with
% parameters CB and SC, instead of conventional Least Squares.
%
% crop specifies border size to crop/ignore  around all sides of the volume. 
%     Should be of the form [ymin xmin zmin; ymax xmax zmax]
%     Default crops 2 pixel border: [2 2; (size(im1) - [2 2])].
%

% default values
if ieNotDefined('robustFlag')
  robustFlag = 0;
end
if ieNotDefined('rotFlag')
  rotFlag = 1;
end
if ieNotDefined('crop')
	crop = [2 2; (size(im1) - [2 2])];
end
if ieNotDefined('CB')
  CB = [];
end
if ieNotDefined('SC')
  SC = [];
end

% Compute derivatives
[fx,fy,ft] = computeDerivatives2(im1,im2);

% Adjust crop
origDims = size(im1);
derivDims = size(fx);
diffDims = origDims - derivDims;
crop(2,:) = crop(2,:) - diffDims;

% Meshgrid on original volume. Then below we throw the edges of the
% meshgrid to make it the same size as the derivative volumes.
% *** Note: This assumes 5tap derivative filters.
[xgrid,ygrid] = meshgrid(1:size(im1,2),1:size(im1,1));
xgrid = xgrid([3:origDims(1)-2],[3:origDims(2)-2]);
ygrid = ygrid([3:origDims(1)-2],[3:origDims(2)-2]);

% Subsample and crop
indicesY = [crop(1,1):2:crop(2,1)];
indicesX = [crop(1,2):2:crop(2,2)];
fx = fx(indicesY,indicesX);
fy = fy(indicesY,indicesX);
ft = ft(indicesY,indicesX);
xgrid = xgrid(indicesY,indicesX);
ygrid = ygrid(indicesY,indicesX);

dimsS=size(fx);
pts=find(~isnan(fx));
pts=find((~isnan(fx))&(~isnan(fy))&(~isnan(ft)));
fx = fx(pts);
fy = fy(pts);
ft = ft(pts);
xgrid = xgrid(pts);
ygrid = ygrid(pts);

if rotFlag
	A= [ fx(:), fy(:), xgrid(:).*fy(:)-ygrid(:).*fx(:)];
else
	A= [xgrid(:).*fx(:), ygrid(:).*fx(:), fx(:),...
		xgrid(:).*fy(:), ygrid(:).*fy(:), fy(:)];
end
b = -ft(:);

if robustFlag
	[p w] = robustMest(A,b,CB,SC);
	w2d = zeros(dimsS);
	w2d(pts)=w;	
else
	p = A\b;
	w2d = [];
end

if rotFlag
	M= [cos(p(3))  -sin(p(3)) p(1);
    	    sin(p(3))  cos(p(3))  p(2);
    	    0     0    1];
else
	M= [1+p(1) p(2)   p(3);
    	    p(4)   1+p(5) p(6);
    	    0      0      1];
end

return;

%%%%%%%%%
% Debug %
%%%%%%%%%

% test with translation
dims=[64 64];
im1=rand(dims);
im2=circularShift(im1,1,0);
% default - rot and LS
estMotion2(im1,im2)
% rot and robust
estMotion2(im1,im2,1,1)
% affine and LS
estMotion2(im1,im2,0,0)
% affine and robust
estMotion2(im1,im2,0,1)

dims=[64 64];
im1=rand(dims);
A= [1 0 .5;
    0 1 .5;
    0 0 1];
im2=warpAffine2(im1,A);
% default - rot and LS
estMotion2(im1,im2)
% rot and robust
estMotion2(im1,im2,1,1)
% affine and LS
estMotion2(im1,im2,0,0)
% affine and robust
estMotion2(im1,im2,0,1)

% test with rotation
dims=[64 64];
in=rand(dims);
theta=atan2(1,max(dims));
A1=[cos(theta) sin(theta) 0;
    -sin(theta) cos(theta) 0;
    0 0 1];
A2=inv(A1);
im1=warpAffine2(in,A1);
im2=warpAffine2(in,A2);
A1*A1
% default - rot and LS
estMotion2(im1,im2)
% rot and robust
estMotion2(im1,im2,1,1)
% affine and LS
estMotion2(im1,im2,0,0)
% affine and robust
estMotion2(im1,im2,0,1)


dims=[64 64];
im1=rand(dims);
A= [1 0 .5;
    0 1 .5;
    0 0 1];
im2=warpAffine2(im1,A);
% default - rot and LS
estMotion2(im1,im2)
% rot and robust
estMotion2(im1,im2,1,1)
% affine and LS
estMotion2(im1,im2,0,0)
% affine and robust
estMotion2(im1,im2,0,1)

%%%%%%%%%%%%%%%%%%%%
% test with outliers
%%%%%%%%%%%%%%%%%%%%

% translation
dims=[64 64];
im1=rand(dims);
A= [1 0 .5;
    0 1 .5;
    0 0 1];
im2=warpAffine2(im1,A);
% putting inconsistent information in upper left corner of im2
Nc=3;
im2(1:round(dims(1)/Nc), 1:round(dims(2)/Nc)) = rand(round(dims(1)/Nc), round(dims(2)/Nc));
% default - rot and LS
ArotLS = estMotion2(im1,im2);
% rot and robust
[ArotRob wr] = estMotion2(im1,im2,1,1);
% affine and LS
AaffLS = estMotion2(im1,im2,0,0);
% affine and robust
[AaffRob wa] = estMotion2(im1,im2,0,1);
ArotLS
ArotRob
AaffLS
AaffRob
imagesc([wr wa]); colormap(gray);axis('image');axis('off')

% test with rotation
dims=[64 64];
in=rand(dims);
theta=atan2(1,max(dims));
A1=[cos(theta) sin(theta) 0;
    -sin(theta) cos(theta) 0;
    0 0 1];
A2=inv(A1);
im1=warpAffine2(in,A1);
im2=warpAffine2(in,A2);
% putting inconsistent information in upper left corner of im2
Nc=3;
im2(1:round(dims(1)/Nc), 1:round(dims(2)/Nc)) = rand(round(dims(1)/Nc), round(dims(2)/Nc));
% default - rot and LS
ArotLS = estMotion2(im1,im2);
% rot and robust
[ArotRob wr] = estMotion2(im1,im2,1,1);
% affine and LS
AaffLS = estMotion2(im1,im2,0,0);
% affine and robust
[AaffRob wa] = estMotion2(im1,im2,0,1);
A=A1*A1
ArotLS
ArotRob
AaffLS
AaffRob
imagesc([wr wa]); colormap(gray);axis('image');axis('off')




