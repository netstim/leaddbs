% getROICoordinates.m
%
%        $Id$
%      usage: scanCoords = getROICoordinates(view,roiNum,<scanNum>,<groupNum>,<baseNum=3>,<straightXform=0>)
%         by: david heeger and justin gardner
%       date: 04/02/07
%    purpose: get roi coordinates in scan coordinates
%             if scanNum is 0, then will compute in the current base
%             coordinates, unless basenum is specified, in which case. 
%             if roinum is a structure, works on the structure
%             rather than the roinum
%             if roinum is a string, will load the roi from
%             the directory
% 
%             If straightXform is set to 1 then the roi is just xform'd
%             without using xfromROIcoords - i.e. it does a straight
%             xform of the coordinates.
function scanCoords = getROICoordinates(view,roiNum,scanNum,groupNum,varargin)

scanCoords = [];
% check arguments
if nargin < 2
  help getROICoordinates
  return
end
straightXform=[];
getArgs(varargin,{'straightXform=0','baseNum=[]'});

% get  scan
if ieNotDefined('scanNum') 
  if ieNotDefined('baseNum')
    scanNum = viewGet(view,'currentScan');
  else
    scanNum = 0;
  end
end

% if roiNum is a string see if it is loaded, otherwise
% try to load it
nROIs = viewGet(view,'numROIs');
if isstr(roiNum)
  if isempty(viewGet(view,'roiNum',roiNum))
    % roi is not loaded in. get it
    view = loadROI(view,roiNum);
    roiNum = viewGet(view,'nROIs');
    if roiNum == nROIs
      disp(sprintf('(getROICoordinates) Could not load ROI'));
      return
    end
  else
    % roi is already installed, just get it
    roiNum = viewGet(view,'roiNum',roiNum);
  end
% if it is a struct use newRoi to set it into the view
elseif isstruct(roiNum)
  [tf roiNum] = isroi(roiNum); 
  if ~tf
    disp(sprintf('(getROICoordinates) Invalid ROI passed in'));
    return
  end
  currentROIname = viewGet(view,'roiName');
  view = viewSet(view,'newROI',roiNum,1);
  roiNum = viewGet(view,'ROINum',roiNum.name);
  %view = viewSet(view,'currentROI',viewGet(view,'ROINum',currentROIname)); %JB: do not change current roi in the view
end

% get the roi transforms
roiVoxelSize = viewGet(view,'roiVoxelSize',roiNum);
roiCoords = viewGet(view,'roiCoords',roiNum);
if isempty(roiCoords),return,end

% make sure we have normalized coordinates
if (size(roiCoords,1)==3)
  roiCoords(4,:) = 1;
end

% get the scan transforms
if scanNum
  if ieNotDefined('groupNum')
    groupNum = viewGet(view,'currentGroup');
  end
  scan2roi = viewGet(view,'scan2roi',roiNum,scanNum,groupNum);
  scanVoxelSize = viewGet(view,'scanVoxelSize',scanNum,groupNum);
else
  % use base xform if scanNum == 0
  %view = viewSet(view,'curGroup',groupNum); %JB: I don't think that's necessary
  if ieNotDefined('baseNum')
    baseNum = viewGet(view,'curbase');
  end
  scan2roi = viewGet(view,'base2roi',roiNum,baseNum);
  scanVoxelSize = viewGet(view,'baseVoxelSize',baseNum);
end  

if (isempty(scan2roi)) 
  disp(sprintf('(getRoiCoordinates) No xform available'));
  return
end

% Use xformROI to supersample the coordinates
if ~straightXform
  scanCoords = round(xformROIcoords(roiCoords,inv(scan2roi),roiVoxelSize,scanVoxelSize));
else
  scanCoords = round(inv(scan2roi)*roiCoords);
end

% return the unique ones
if ~isempty(scanCoords)
  scanCoords = unique(scanCoords','rows')';
  scanCoords = scanCoords(1:3,:);
end

if ~isempty(scanCoords) && (scanNum ~= 0)

  % check scan dimensions
  scanDims = viewGet(view,'dims',scanNum,groupNum);

  % make sure we are inside scan dimensions
  xCheck = (scanCoords(1,:) >= 1) & (scanCoords(1,:) <= scanDims(1));
  yCheck = (scanCoords(2,:) >= 1) & (scanCoords(2,:) <= scanDims(2));
  sCheck = (scanCoords(3,:) >= 1) & (scanCoords(3,:) <= scanDims(3));

  % only return ones that are in bounds
  scanCoords = scanCoords(:,find(xCheck & yCheck & sCheck));
end



