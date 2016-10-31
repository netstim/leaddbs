function thisView = restrictROI(thisView,ROInum,scan)

% function thisView = restrictROI(thisView,[ROInum],[scan])
%
% ROInum can be either a name, a number, or empty (for current ROI).
% Default: current ROI
%
% scan can be either a number or empty (for current scan)
% Default: current scan
%
% djh 9/2005

if ieNotDefined('ROInum')
  ROInum = viewGet(thisView,'currentroinum');
end
if ischar(ROInum)
  ROInum = viewGet(thisView,'roinum',ROInum);
end

if ieNotDefined('scan')
  scan = viewGet(thisView,'curscan');
end

ROIcoords = viewGet(thisView,'roiCoords',ROInum);
% Save prevCoords for undo
thisView = viewSet(thisView,'prevROIcoords',ROIcoords);

% Transform ROI roiScanCoords to overlay
roiScanCoords = round( viewGet(thisView,'scan2roi',ROInum,scan) \ ROIcoords); 

coordsInfo.base2overlay = eye(4);
coordsInfo.baseCoordsHomogeneous = roiScanCoords;
coordsInfo.baseDims = [size(ROIcoords,2) 1 1];

%find which voxels are not clipped in the current overlay(s) and overlayAlpha (in overlay space)
overlayList  = viewGet(thisView,'curOverlay');
nOverlays = length(overlayList);
cOverlay=0;
alphaOverlayList = zeros(size(overlayList));
for iOverlay=overlayList
  cOverlay = cOverlay+1;
  thisAlphaOverlay = viewGet(thisView,'overlayNum',viewGet(thisView,'alphaOverlay',iOverlay));
  if ~isempty(thisAlphaOverlay)
    alphaOverlayList(cOverlay) = thisAlphaOverlay;
  end
end
roiMask = maskOverlay(thisView,[overlayList alphaOverlayList],scan,coordsInfo);
roiMask = reshape(roiMask{1},[size(roiMask{1},1) nOverlays, 2]);
% keep the corresponding voxels in ROI space
cOverlay=0;
for iOverlay=overlayList
  cOverlay = cOverlay+1;
  if ~alphaOverlayList(cOverlay) %if there is no alpha overlay, it's as if all voxels were 1
    roiMask(:,cOverlay,2)=1;
  end
end    
%Keep voxels that are non-zero in any of the overlays, but non-zero both in overlay and alphaOverlay, 
ROIcoords = ROIcoords(:,any(all(roiMask,3),2));


thisView = viewSet(thisView,'roiCoords',ROIcoords,ROInum);

return

%old code, now replaced by call to maskOverlay
for overlayNum = 1:viewGet(thisView,'numberOfOverlays')
  overlayData = viewGet(thisView,'overlayData',scan,overlayNum);
  overlayClip = viewGet(thisView,'overlayClip',overlayNum);
  if ~isempty(overlayData) & ~isempty(ROIcoords)
    % Transform ROI roiScanCoords to overlay
    xform = inv(viewGet(thisView,'scan2roi',ROInum,scan));
    roiScanCoords = round(xform * ROIcoords); 
    overlaySize = size(overlayData);
    % Find roiScanCoords that are within the overlay size
    for v = 1:3
      indices = find((roiScanCoords(v,:) > 0) & (roiScanCoords(v,:) <= overlaySize(v)));
      roiScanCoords = roiScanCoords(:,indices);
      ROIcoords = ROIcoords(:,indices);
    end
    % Find overlay voxels that are within clip
    indices = sub2ind(overlaySize,roiScanCoords(1,:),roiScanCoords(2,:),roiScanCoords(3,:));
    data = overlayData(indices);
    if diff(overlayClip) > 0
      indices = find(data >= overlayClip(1) & data <= overlayClip(2));
    else
      indices = find(data >= overlayClip(1) | data <= overlayClip(2));
    end
    ROIcoords = ROIcoords(:,indices);
  else
    ROIcoords = [];
  end
  thisView = viewSet(thisView,'roiCoords',ROIcoords,ROInum);
end

return;
