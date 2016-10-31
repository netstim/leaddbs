% mlrImageOrient.m
%
%      usage: [data h xform] = mlrImageOrient(orient,data,h)
%         by: justin gardner
%       date: 09/05/11
%    purpose: re-orient image to a standard orientation like LPI
%             This will fix the header to have the correct
%             qform/sform
%
%       e.g.: [d h] = mlrImageOrient('LPI',data,h);
%        
%             you can also get the xform that was applied to reorient the image
%             [d h xform] = mlrImageOrient('LPI',data,h);
%  
%
function [data h xform] = mlrImageOrient(orient,varargin)

% check arguments
if nargin < 2
  help mlrImageOrient
  return
end

[imageArgs] = mlrImageParseArgs(varargin);
if isempty(imageArgs) || ~mlrImageIsImage(imageArgs{1})
  disp(sprintf('(mlrImageOrient) No image to reorient'));
  return
else
  [data h] = mlrImageLoad(imageArgs{1});
end

% convert the orientation string into something easier to handle
[tf axisIndex] = ismember(orient,'LRPAIS');
desiredAxis = ceil(axisIndex/2);
desiredDir = 2*isodd(axisIndex)-1;
[dummy desiredAxisReverse] = sort(desiredAxis);

% check to see if we have all axis represented
if length(unique(desiredAxis)) ~= 3
  disp(sprintf('(mlrImageOrient) Invalid orientation, must have one each of the 3 axis: L/R P/A I/S'));
  return
end

% check header
if ~mlrImageIsHeader(h)
  disp(sprintf('(mlrImageOrient) Header is not a standard mlrImage header'));
  return
end

% check for qform
if isempty(h.qform)
  disp(sprintf('(mlrImageOrient) No qform available in image header'));
  return
end

% get the axis directions
axisLabels = mlrImageGetAxisLabels(h.qform);

% now make axisMapping and axisReverseMapping appropriate
% for the desired orientation
axisMapping = desiredAxisReverse(axisLabels.mapping);
axisReverseMapping = axisLabels.reverseMapping(desiredAxis);
desiredDir = desiredDir(axisMapping);
axisDir = axisLabels.dir;

% flip each axis that goes in the opposite direction
for i = 1:3
  if axisDir(i)*desiredDir(i) == -1
    data = flipdim(data,i);
  end
end

% get the permutation order
permutationOrder = 1:h.nDim;
permutationOrder(1:3) = axisReverseMapping;

% and permute
data = permute(data,permutationOrder);

% make the xform matrix
for i = 1:3
  xform(i,1:4) = 0;
  xform(i,axisMapping(i)) = axisDir(i)*desiredDir(i);
  if axisDir(i)*desiredDir(i) == -1
    xform(i,4) = h.dim(i)-1;
  end
end
xform(4,:) = [0 0 0 1];

% fix qform and sform
h.qform = h.qform*xform;
if ~isempty(h.sform)
  h.sform = h.sform*xform;
end

% fix the pixdim
if length(h.pixdim) >= 3
  h.pixdim(1:3) = abs((xform(1:3,1:3)*h.pixdim(1:3)')');
end

% fix dimensions
h.dim = size(data);