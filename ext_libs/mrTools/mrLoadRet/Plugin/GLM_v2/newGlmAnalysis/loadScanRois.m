% loadScanRois - load timeseries corresponding to the voxels of one or several loaded Rois
%
%        $Id$
%      usage: [d, roiVoxelIndices, roiCoords] = loadScanRois(thisView,<scanList>,<roiList>)
%         by: julien besle
%       date: 2010-04-26
%     inputs: thisView
%             scanList: list of scans whose times-series will be concatenated (default = current scan)
%             roiList: list of ROIs to be loaded (default = current ROI(s))
%    outputs: d is a data structure with all the voxels on the first dimensions and time-series on the 4th dimension
%             roiVoxelIndices is a cell array indexing which voxel in the data are part of each roi
%             roiCoords is a cell array with the voxel coordinates in the scan
%             If no scan/roi are specified, currently loaded scan/roi(s) are used
%
%    purpose: load timeseries corresponding to the voxels of one or several Rois (wrapper around loadScan.m)
%

function [d, roiVoxelIndices, roiCoords] = loadScanRois(thisView,scanList,roiList)

if ieNotDefined('scanList'),scanList = viewGet(thisView,'curScan');end
if ieNotDefined('roiList'),roiList = viewGet(thisView,'currentROI');end

%get the smallest box in the scan including all the rois in the list
[subsetBox, whichRoi]  = getRoisBox(thisView,scanList,0,roiList);  
isInARoi = find(any(whichRoi,4));
% isInARoi = any(whichRoi,4);

precision = mrGetPref('defaultPrecision');
cScan=0;
for iScan=scanList
  cScan = cScan+1;
  dummy = loadScan(thisView,iScan,[],subsetBox(3,:),precision,subsetBox(1,:),subsetBox(2,:));
  %only keep the data that's in at least one of the ROIs
  dummy.data = permute(dummy.data,[4 1 2 3]);
  dummy.data = dummy.data(:,isInARoi);
  if cScan==1
    d = dummy;
  else
    d.data = cat(1,d.data,dummy.data);
  end
end
clear('dummy');

%checking for NaNs in the tseries
isnanVoxels = any(isnan(d.data),1);
if nnz(isnanVoxels)
   mrWarnDlg('(loadScanRois) removing undefined time-series');
   d.data = d.data(:,~isnanVoxels);
   isnanIndex = isInARoi(isnanVoxels);
   isInARoi = isInARoi(~isnanVoxels);
else
  isnanIndex=[];
end
d.data = permute(d.data,[2 3 4 1]);
d.dim = size(d.data);

if nargout>1
  roiVoxelIndices=cell(1,length(roiList));
end
if nargout>2
  roiCoords=cell(1,length(roiList));
end

%for each roi, identify undefined voxels and find the indices of the voxels in the loaded data
subsetDims = (diff(subsetBox,1,2)+1)';
for iRoi = 1:length(roiList)
  isInThisRoi = find(whichRoi(:,:,:,iRoi));
  nansInThisRoi = intersect(isInThisRoi,isnanIndex);
  if ~isempty(nansInThisRoi)
    fprintf(1,['ROI ' num2str(iRoi) ' (' viewGet(thisView,'roiName',roiList(iRoi)) '):\n']);
    [isnanXCoords,isnanYCoords,isnanZCoords] = ind2sub(subsetDims,nansInThisRoi);
    disp((repmat(subsetBox(:,1)'-1,size(isnanXCoords,1),1)+[isnanXCoords,isnanYCoords,isnanZCoords])');
    isInThisRoi = setdiff(isInThisRoi,nansInThisRoi);
  end
  if nargout>1
    [boxIndices, roiVoxelIndices{iRoi}] = intersect(isInARoi,isInThisRoi);
  end
  if nargout>2
    [x y z] = ind2sub(subsetDims,boxIndices);
    roiCoords{iRoi} = zeros(3,length(boxIndices));
    roiCoords{iRoi}(1,:)=x+subsetBox(1,1)-1;
    roiCoords{iRoi}(2,:)=y+subsetBox(2,1)-1;
    roiCoords{iRoi}(3,:)=z+subsetBox(3,1)-1;
  end
end
