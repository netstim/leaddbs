function pasteOverlay(thisView, overlays)
% pasteOverlay.m
%
%        $Id$
%      usage: pasteOverlayBaseCoords(thisView, clipboard)
%         by: julien besle
%       date: 25/01/2010
%    purpose: pastes one or many MLR overlay(s) from the clipboard, and converts them to the current scan space

if ~isanalysis(viewGet(thisView,'analysis'))
    mrWarnDlg('(pasteOverlay) Overlays must be pasted into an analysis. Use Edit -> Analysis -> New Analysis.')
    return
end

nScans = viewGet(thisView,'nScans');
scanList = 1:nScans;
if isempty(scanList)
   mrWarnDlg('(pasteOverlay) Could not paste overlay because the group is empty (no scan)' );
   return;
end
while length(scanList)~= length(overlays(1).data)
   scanList=selectInList(thisView,'scans');
   if isempty(scanList)
     return;
   end
end
curGroupName = viewGet(thisView,'groupName');
curGroupNum =  viewGet(thisView,'groupNum',curGroupName);

for iOverlay = 1:length(overlays)
   [check thisOverlay] = isoverlay(overlays(iOverlay));
   if ~check
       mrErrorDlg('(pasteOverlay) Cannot paste. Clipboard does not contain a valid overlay. Use Edit -> Overlay -> Copy Overlay.')
   end
   fromGroupNum = viewGet(thisView,'groupNum',thisOverlay.groupName);
   scan2scan = viewGet(thisView,'scan2scan',1,curGroupNum,1,fromGroupNum);
   scandims = viewGet(thisView, 'scandims', 1);
   overlays(iOverlay).data = cell(1,nScans);
   if ~all(all(abs(scan2scan - eye(4))<1e-6)) %check if we're in the scan space
      %transform values in current scan space
      [Ycoords,Xcoords,Zcoords] = meshgrid(1:scandims(2),1:scandims(1),1:scandims(3));
      for iScan = 1:length(scanList)
         fprintf(1,['Overlay ' num2str(iOverlay) ', Scan ' num2str(iScan) '\n']);
         overlays(iOverlay).data{scanList(iScan)} = getNewSpaceOverlay(thisOverlay.data{iScan}, scan2scan,Xcoords,Ycoords,Zcoords);
      end
   else
      overlays(iOverlay).data(scanList) = thisOverlay.data;
   end
   overlays(iOverlay).groupName =  curGroupName;
   overlays(iOverlay).name = [overlays(iOverlay).name ' (from ' thisOverlay.groupName ')'];
   if isfield(overlays(iOverlay),'alphaOverlay') && ~isempty(overlays(iOverlay).alphaOverlay)
     overlays(iOverlay).alphaOverlay = [overlays(iOverlay).alphaOverlay ' (from ' thisOverlay.groupName ')'];
   end
   %remove any specific reconcile function and params just in case, this overlay has to emancipate 
   overlays(iOverlay).reconcileFunction = 'defaultReconcileParams';
   overlays(iOverlay).params = [];

end
thisView = viewSet(thisView,'newOverlay',overlays);


