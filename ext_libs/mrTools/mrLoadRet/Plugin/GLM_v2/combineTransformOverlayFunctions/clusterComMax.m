
function [clusterOverlay, comOverlay, maxOverlay] = clusterComMax(inputOverlay, test, minSize)
% [clusterOverlay] = clusterMax(inputOverlay, test, minSize)
%
%   finds clusters of connected voxels and their max and center of Mass
%
% Input:
%     overlay 3D data array
%
% One Mandatory additional scalar input
%     test:  anonymous function to threshold the map (e.g. @(x)x>1.65)
%     minSize (optional): minimum number of voxel in cluster (default=0)
%        
% Output: select 3 output overlays
%     - clusterOverlay 3D array where cluster values represented by the number of scan voxels in the cluster
%     - center of mass overlays, were only the COM of the cluster are left
%     - max overlay, with only the max value of each cluster
%
% jb 21/03/2011
%
% $Id$ 

if ~ismember(nargin,[2 3])
  help clusterMax;
  return;
elseif ~isa(test,'function_handle')
   mrWarnDlg('(clusterMax) first additional argument must be an anonymous function');
   return;
end

if ieNotDefined('minSize')
  minSize=0;
end


clusterOverlay = test(inputOverlay);

%find contiguous voxels
cc = bwconncomp(clusterOverlay,6);
sizeObjects = zeros(1,cc.NumObjects);
for iObject = 1:cc.NumObjects
  sizeObjects(iObject) = length(cc.PixelIdxList{iObject});
end
minSizeObjects = find(sizeObjects>minSize);
sizeObjects = sizeObjects(minSizeObjects);
cc.NumObjects = length(sizeObjects);
cc.PixelIdxList = cc.PixelIdxList(minSizeObjects);
[dump,orderedSizes] = sort(sizeObjects,'descend');
cc.PixelIdxList = cc.PixelIdxList(orderedSizes);
sizeObjects = sizeObjects(orderedSizes);

clusterOverlay = NaN(size(clusterOverlay));
comOverlay = NaN(size(clusterOverlay));
maxOverlay = NaN(size(clusterOverlay));
for iObject = 1:cc.NumObjects
  clusterOverlay(cc.PixelIdxList{iObject})=sizeObjects(iObject);
  clusterData = inputOverlay(cc.PixelIdxList{iObject});
  
  %find the center of mass
  [coordsX,coordsY,coordsZ] = ind2sub(size(clusterOverlay),cc.PixelIdxList{iObject});
  COM = round(sum(repmat(clusterData,1,3)/sum(clusterData) .* [coordsX,coordsY,coordsZ],1));
  comOverlay(COM(1),COM(2),COM(3))=sizeObjects(iObject);
  
  %find the max 
  [maxData,maxIndex] = max(clusterData);
  maxOverlay(cc.PixelIdxList{iObject}(maxIndex))=maxData;
  
end

