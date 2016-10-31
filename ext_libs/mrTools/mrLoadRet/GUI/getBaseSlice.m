% getBaseSlice.m
%
%        $Id$ 
%      usage: [baseIm baseCoords baseCoordsHomogeneous] = getBaseSlice(view, sliceNum, sliceIndex, rotate, baseNum, baseType)
%         by: justin gardner
%       date: 03/05/10
%    purpose: DJH - pulled out from refreshMLRDisplay for use in other functions -jg
%             extracts base image and corresponding coordinates
function [baseIm,baseCoords,baseCoordsHomogeneous] = getBaseSlice(view,sliceNum,sliceIndex,rotate,baseNum,baseType)

baseCoords = [];
baseCoordsHomogeneous = [];
baseIm = [];

% viewGet
volSize = viewGet(view,'baseDims',baseNum);
baseData = viewGet(view,'baseData',baseNum);

% get the crop type
% for regular images we want loose, for
% flat maps, we want crop
if baseType == 0,cropType = 'loose';else cropType = 'crop';end

if ~isempty(volSize)

  % Generate coordinates with meshgrid
  switch sliceIndex
     case 1
      % validate sliceNum against base dimensions
      if sliceNum>volSize(1),sliceNum = volSize(1);end  
      % set up the coordinates
      x = sliceNum * ones(volSize(2),volSize(3));
      [z,y] = meshgrid(1:volSize(3),1:volSize(2));
    case 2
      if sliceNum>volSize(2),sliceNum = volSize(2);end  
      y = sliceNum * ones(volSize(1),volSize(3));
      [z,x] = meshgrid(1:volSize(3),1:volSize(1));
    case 3
      if sliceNum>volSize(3),sliceNum = volSize(3);end  
      z = sliceNum * ones(volSize(1),volSize(2));
      [y,x] = meshgrid(1:volSize(2),1:volSize(1));
  end

  % if there are baseCoords associated with this base
  % then use those instead of image coordinates (this
  % is the case for flat files)
  baseCoordMap = viewGet(view,'baseCoordMap',baseNum);
  if ~isempty(baseCoordMap)
    x = permute(baseCoordMap.coords(:,:,1,1,:),[1 2 5 3 4]);
    y = permute(baseCoordMap.coords(:,:,1,2,:),[1 2 5 3 4]);
    z = permute(baseCoordMap.coords(:,:,1,3,:),[1 2 5 3 4]);
  end

  % Rotate coordinates
  if (rotate ~= 0)
    for iDepth = 1:size(x,3) %this is if we're taking several depth bins in a flat map (there will only be one "depth" for volumes)
      % this is nearest neighbor interpolation for a reason - if we draw ROIs on this, we want
      % the coordinates to match no matter what rotation is specified. If we use nearest, then we
      % always keep actualy coordinates from the base anatomy. If we interpolate, all bets are off
      % and you get little broken pixels in the ROI when you rotate (or convertROICorticalDepth)
      newX(:,:,iDepth) = imrotate(x(:,:,iDepth),rotate,'nearest',cropType);
      newY(:,:,iDepth) = imrotate(y(:,:,iDepth),rotate,'nearest',cropType);
      newZ(:,:,iDepth) = imrotate(z(:,:,iDepth),rotate,'nearest',cropType);
    end
    x=newX;
    y=newY;
    z=newZ;
  end

  % Reformat base coordinates
  imageDims = [size(x,1) size(x,2)];
  numDepths = size(x,3);
  numPixels = prod(imageDims);
  xvec = reshape(x,1,numPixels*numDepths);
  yvec = reshape(y,1,numPixels*numDepths);
  zvec = reshape(z,1,numPixels*numDepths);
  baseCoordsHomogeneous = [xvec; yvec; zvec; ones(1,numPixels*numDepths)];
  baseCoordsHomogeneous = reshape(baseCoordsHomogeneous,[4 numPixels numDepths]);
  baseCoords = permute(reshape(baseCoordsHomogeneous(1:3,:),[3 imageDims numDepths]),[2 3 1 4]);
end

% Extract base image
if ~isempty(baseData)
  switch sliceIndex
    case 1
      baseIm = squeeze(baseData(sliceNum,:,:));
    case 2
      baseIm = squeeze(baseData(:,sliceNum,:));
    case 3
      baseIm = squeeze(baseData(:,:,sliceNum));
  end
  if (rotate ~= 0)
    baseIm = imrotate(baseIm,rotate,'bilinear',cropType);
  end
end

%disp(sprintf('(getBaseSlice) DEBUG: sliceIndex=%i, sliceNum=%i',sliceIndex,sliceNum));
return
