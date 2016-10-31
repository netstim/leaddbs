% mlrImageXform.m
%
%      usage: [d h] = mlrImageXform(image,commands)
%         by: justin gardner
%       date: 09/10/11
%    purpose: transforms an image while also adjusting header info including
%             xformation matrices accordingly. command is any of
%             the following:
%
%             swapXY,swapXZ,swapYZ
%             flipX,flipY,flipZ
%
%             shiftX=n,shiftY=n,shiftZ=n
%             rotateXY=deg,rotateYZ=deg,rotateXZ=deg
%
%       e.g.: [d h] = mlrImageXform('image.hdr','swapXY=1','flipZ=1','shiftX=13.2',rotateXY=30');
%
%
function [d h xform] = mlrImageXform(varargin)

% check arguments
if nargin < 1
  help mlrImageXform
  return
end
d = [];h= [];xform=[];

% parse arguments
[imageArgs otherArgs] = mlrImageParseArgs(varargin);
verbose=[];
swapXY=[];swapXZ=[];swapYZ=[];
flipX=[];flipY=[];flipZ=[];
shiftX=[];shiftY=[];shiftZ=[];
rotateXY=[];rotateXZ=[];rotateYZ=[];
xMin=[];xMax=[];yMin=[];yMax=[];zMin=[];zMax=[];
interpMethod=[];
applyToHeader=[];applyToData=[];
getArgs(otherArgs,{'verbose=1','swapXY=0','swapXZ=0','swapYZ=0','flipX=0','flipY=0','flipZ=0','shiftX=0','shiftY=0','shiftZ=0','rotateXY=0','rotateXZ=0','rotateYZ=0','interpMethod=[]','applyToHeader=1','applyToData=1','xMin=1','xMax=inf','yMin=1','yMax=inf','zMin=1','zMax=inf'});

% check that we have an image to xform
if length(imageArgs) < 1
  disp(sprintf('(mlrImageXform) Must specify atleast one image to xform'));
  return
end

% get default interpMethod for xforms that use interp3
if isempty(interpMethod)
  interpMethod = mrGetPref('interpMethod');
end

% init xform
xform = eye(4);

% now cycle through images xforming them and putting in output
for iImage = 1:length(imageArgs)
  % load the image
  [d h] = mlrImageLoad(imageArgs{iImage},'returnRaw=1');
  if isempty(d)
    disp(sprintf('(mlrImageXform) Could not open image: %s',mlrImageArgFilename(imageArgs{iImage})));
  else
    % apply swapXY
    if swapXY
      if verbose,disp(sprintf('(mlrImageXform) Swapping XY'));end
      % swap the data
      if applyToData
	permuteDims = 1:h.nDim;
	permuteDims(1:2) = [2 1];
	d = permute(d,permuteDims);
      end
      % swap the header
      [h xform] = applyXform([0 1 0 0;1 0 0 0;0 0 1 0;0 0 0 1],h,d,applyToHeader,xform);
    end
    % apply swapXZ
    if swapXZ
      if verbose,disp(sprintf('(mlrImageXform) Swapping XZ'));end
      % swap the data
      if applyToData
	permuteDims = 1:h.nDim;
	permuteDims(1:3) = [3 2 1];
	d = permute(d,permuteDims);
      end
      % swap the header
      [h xform] = applyXform([0 0 1 0;0 1 0 0;1 0 0 0;0 0 0 1],h,d,applyToHeader,xform);
    end
    % apply swapYZ
    if swapYZ
      if verbose,disp(sprintf('(mlrImageXform) Swapping YZ'));end
      % swap the data
      if applyToData
	permuteDims = 1:h.nDim;
	permuteDims(2:3) = [3 2];
	d = permute(d,permuteDims);
      end
      % apply to header
      [h xform] = applyXform([1 0 0 0;0 0 1 0;0 1 0 0;0 0 0 1],h,d,applyToHeader,xform);
    end
    % apply flipX
    if flipX
      if verbose,disp(sprintf('(mlrImageXform) Flipping X'));end
      % flip data
      if applyToData
	d = flipdim(d,1);
      end
      % apply to header
      [h xform] = applyXform([-1 0 0 h.dim(1)-1;0 1 0 0;0 0 1 0;0 0 0 1],h,d,applyToHeader,xform);
    end
    % apply flipY
    if flipY
      if verbose,disp(sprintf('(mlrImageXform) Flipping Y'));end
      % flip data
      if applyToData
	d = flipdim(d,2);
      end
      % apply to header
      [h xform] = applyXform([1 0 0 0;0 -1 0 h.dim(2)-1;0 0 1 0;0 0 0 1],h,d,applyToHeader,xform);
    end
    % apply flipZ
    if flipZ
      if verbose,disp(sprintf('(mlrImageXform) Flipping Z'));end
      % flip data
      if applyToData
	d = flipdim(d,3);
      end
      % apply to header
      [h xform] = applyXform([1 0 0 0;0 1 0 0;0 0 -1 h.dim(3)-1;0 0 0 1],h,d,applyToHeader,xform);
    end
    % apply shifts
    if any([shiftX shiftY shiftZ])
      if verbose,disp(sprintf('(mlrImageXform) Using %s to shift by: [%s]',interpMethod,mlrnum2str([shiftX shiftY shiftZ])));end
      % flip data
      if applyToData
	x = 1+shiftX:h.dim(1)+shiftX;
	y = 1+shiftY:h.dim(2)+shiftY;
	z = 1+shiftZ:h.dim(3)+shiftZ;
	[x y z] = ndgrid(x,y,z);
	d = interpn(d,x,y,z,interpMethod);
      end
      % apply to header
      [h xform] = applyXform([1 0 0 shiftX;0 1 0 shiftY;0 0 1 shiftZ;0 0 0 1],h,d,applyToHeader,xform);
    end
    % apply rotation
    if any([rotateXY rotateXZ rotateYZ])
      if verbose,disp(sprintf('(mlrImageXform) Using %s to rotate: %s',interpMethod,mlrnum2str([rotateXY rotateXZ rotateYZ])));end
      % flip data
      if applyToData
	% shift coordinates to center
	x = (1:h.dim(1))-h.dim(1)/2;
	y = (1:h.dim(2))-h.dim(2)/2;
	z = (1:h.dim(3))-h.dim(3)/2;
	% make into grid, note the swapping of x and y
	[x y z] = ndgrid(x,y,z);
	% make into homogenous coordinates
	coords(1,:) = x(:);
	coords(2,:) = y(:);
	coords(3,:) = z(:);
	coords(4,:) = 1;
	% now multiply against the rotation matrix
 	rotMatrix = makeRotMatrix3D(pi*rotateXZ/180,pi*rotateYZ/180,pi*rotateXY/180,[0 0 0]);
	coords = rotMatrix*coords;
	% and shift back from center
	coords(1,:) = coords(1,:) + h.dim(1)/2;
	coords(2,:) = coords(2,:) + h.dim(2)/2;
	coords(3,:) = coords(3,:) + h.dim(3)/2;
	% and interp
	d = interpn(d,coords(1,:),coords(2,:),coords(3,:),interpMethod);
	d = reshape(d,h.dim(1),h.dim(2),h.dim(3));
      end
      % apply to header
      shiftToCenter = [1 0 0 h.dim(1)/2;0 1 0 h.dim(2)/2;0 0 1 h.dim(3)/2;0 0 0 1];
      [h xform] = applyXform(shiftToCenter*rotMatrix*inv(shiftToCenter),h,d,applyToHeader,xform);
    end
  end
  % do any needed adjustment of dimensions
  [d h xform] = adjustDims(d,h,xMin,xMax,yMin,yMax,zMin,zMax,applyToHeader,applyToData,xform);
  % make sure dim is correct
  h.dim = size(d);
  % fix pixdim - thus if we have anisotropic voxels and we swap
  % dims or rotate this will fix the dimensions to be appropriate
  if applyToHeader
    if length(h.pixdim) >= 3
      h.pixdim(1:3) = abs((xform(1:3,1:3)*h.pixdim(1:3)')');
    end
  end
end

%%%%%%%%%%%%%%%%%%%%
%%   applyXform   %%
%%%%%%%%%%%%%%%%%%%%
function [h totalXform] = applyXform(xform,h,d,applyToHeader,totalXform)

% nothing to do if we are not applying to header
if ~applyToHeader,return,end

% keep track of what the total tranform is
totalXform = totalXform*xform;

% apply to qform if it is present
if ~isempty(h.qform)
  h.qform = h.qform*xform;
end

% apply to sform if it is present
if ~isempty(h.sform)
  h.sform = h.sform*xform;
end

% get the dimensions of the data
sized = size(d);
h.dim(1:length(sized)) = sized;

%%%%%%%%%%%%%%%%%%%%
%%   adjustDims   %%
%%%%%%%%%%%%%%%%%%%%
function [data h totalXform] = adjustDims(data,h,xMin,xMax,yMin,yMax,zMin,zMax,applyToHeader,applyToData,totalXform)

% if nothing to do, just return
if (xMin == 1) && (xMax == inf) && (yMin == 1) && (yMax == inf) && (zMin == 1) && (zMax == inf)
  return
end

% get current dimensions
dims = size(data);

% make the dimensions valid
xMin = round(max(1,xMin));
xMax = round(min(xMax,dims(1)));
yMin = round(max(1,yMin));
yMax = round(min(yMax,dims(2)));
zMin = round(max(1,zMin));
zMax = round(min(zMax,dims(3)));

% make some checks
if (zMax < zMin) || (yMax < yMin) || (xMax < xMin)
  disp(sprintf('(mlrImageXform:adjustDims) Could not adjust dims to [%i:%i,%i:%i,%i:%i]',xMin,xMax,yMin,yMax,zMin,zMax));
  return
end

% adjust data size
if applyToData
  data = data(xMin:xMax,yMin:yMax,zMin:zMax,:);
end

% find out how much the new xMin, yMin, zMin have translated the image
t = h.qform * [xMin yMin zMin 1]' - h.qform * [1 1 1 1]';
t = [zeros(3,3) t(1:3,1); 0 0 0 0];
if applyToHeader
  h.qform = t+h.qform;
end
totalXform = totalXform + t;

% reset the dims
h.dim = size(data);



