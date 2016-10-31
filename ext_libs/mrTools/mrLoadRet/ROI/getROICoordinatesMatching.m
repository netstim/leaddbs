% getROICoordinatesMatching.m
%
%      usage: scanCoords = getROICoordinatesMatching(v,roiNum,scanNum,matchScanNum,<groupNum>,<matchGroupNum>)
%         by: justin gardner
%       date: 01/21/09
%    purpose: This gets the rois scan coordinates for the scanNum/groupNum.
%             It does this in a special way. It first gets the scan coordinates of
%             the roi for the matchScanNum/matchGroupNum. It then finds the corresponding
%             coordinate for each one of these coordinates in the scanNum. This is done
%             so that the roi for the scanNum has exactly the same voxels as the
%             matchScanNum roi. This can be useful for classification protocols in which
%             you want to have the same voxels from different sessions. If matchGroupNum
%             is not specified it is assumed to be the same as the groupNum. This would
%             usually be run on two scans that have been imported from different sessions
function scanCoords = getROICoordinatesMatching(v,roiNum,scanNum,matchScanNum,groupNum,matchGroupNum)

scanCoords = [];
% check arguments
if ~any(nargin == [4 5 6])
  help getROICoordinatesMatching
  return
end

% get default groups
if ieNotDefined('groupNum')
  groupNum = viewGet(v,'currentGroup');
end
if ieNotDefined('matchGroupNum')
  matchGroupNum = groupNum;
end

% get the source coordinates
matchScanCoords = getROICoordinates(v,roiNum,matchScanNum,matchGroupNum);
if size(matchScanCoords,1) == 3
  matchScanCoords(4,:) = 1;
end

% get the scan2scan transforms
scan2scan = viewGet(v,'scan2scan',matchScanNum,matchGroupNum,scanNum,groupNum);

disp(sprintf('(getROICoordinatesMatching) Matching %s:%i ROI coordinates to %s:%i',viewGet(v,'groupName',groupNum),scanNum,viewGet(v,'groupName',matchGroupNum),matchScanNum));

% cycle over voxels
for voxNum = 1:size(matchScanCoords,2);
  scanCoords(:,voxNum) = round(scan2scan*matchScanCoords(:,voxNum));
%  disp(sprintf('(%s) %s:%i [%i %i %i] -> %s:%i [%i %i %i]',mfilename,viewGet(v,'groupName',matchGroupNum),matchScanNum,matchScanCoords(1,voxNum),matchScanCoords(2,voxNum),matchScanCoords(3,voxNum),viewGet(v,'groupName',groupNum),scanNum,scanCoord(1,voxNum),scanCoord(2,voxNum),scanCoord(3,voxNum)));
end


