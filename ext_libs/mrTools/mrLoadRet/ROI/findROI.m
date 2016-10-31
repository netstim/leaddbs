% findROI.m
%
%        $Id$
%      usage: findROI()
%         by: justin gardner
%       date: 09/25/07
%    purpose: looks for an ROI in the base anatomy (is called form mrLoadRetGUI)
%
function retval = findROI(view)

retval = [];

% check arguments
if ~any(nargin == [1])
  help findROI
  return
end

if viewGet(view,'baseType') ~= 0
  disp(sprintf('(findROI) Base anatomy cannot be flat or surface to look for an ROI'));
  return
end
% get the roi
roiNum = viewGet(view,'currentROI');
if isempty(roiNum),return,end

baseNum = viewGet(view,'currentBase');
baseName = fixBadChars(viewGet(view,'baseName',baseNum));
sliceIndex = viewGet(view,'baseSliceIndex',baseNum);

% get roi coordinates from the cache
for r = viewGet(view,'currentROI')
  % look in cache for roi
  roi = viewGet(view,'ROICache',r);

  if ~isempty(roi)
    if ((~isfield(roi,baseName)) || ...
        (length(roi.(baseName)) < sliceIndex) || ...
        (isempty(roi.(baseName){sliceIndex})))
      msgbox(sprintf('ROI %s has no voxels in the current anatomy',viewGet(view,'roiName',roiNum)));
    else
      s = roi.(baseName){sliceIndex}.s;
      break;
    end
  end
end

if isempty(roi)
  disp(sprintf('(findROI) ROI cache is emtpy (probably you are not viewing ROIs)'));
  return
end

% get current slice
curSlice = viewGet(view,'curSlice');

% find the closest slice
distanceToCurrentSlice = abs(s-curSlice);
closestSlice = s(first(find(min(distanceToCurrentSlice)==distanceToCurrentSlice)));

% set the slice
if ~isempty(closestSlice)
  viewSet(view,'curSlice',closestSlice);
  refreshMLRDisplay(view.viewNum);
else
  mrWarnDlg(sprintf('(findROI) Could not find ROI %i: %s in base anatomy',roiNum,viewGet(view,'roiname',roiNum)));
end