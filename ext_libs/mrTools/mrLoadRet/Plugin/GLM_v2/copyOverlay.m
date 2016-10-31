% copyOverlay.m
%
%      usage: overlays = copyOverlayBaseCoords(thisView)
%         by: julien besle, based on mrExprt2SR by eli merriam
%       date: 03/20/07
%    purpose: copies many MLR overlays to the clipboard
%        $Id$

function overlays = copyOverlay(thisView)

keepAsking = 1;
overlays = [];

while keepAsking
   scanList = 1:viewGet(thisView,'nScans');
   if length(scanList)>1
      scanList = selectInList(thisView,'scans');
   end
   if isempty(scanList)
      return;
   else
      keepAsking = 0;

      overlayList = 1:viewGet(thisView,'nOverlays');
      if length(overlayList)>1
         overlayList = selectInList(thisView,'overlays');
         if isempty(overlayList)
            if viewGet(thisView,'nScans')==1
               return
            else
               keepAsking=1;
            end
         end
      end

   end
end

overlayNames = viewGet(thisView,'overlayNames');
overlayNames = overlayNames(overlayList);
cOverlay = 0;
for iOverlay = 1:length(overlayList)
  cOverlay = cOverlay+1;
  overlays = copyFields(viewGet(thisView,'overlay',overlayList(iOverlay)),overlays,cOverlay);
  overlays(cOverlay).data = overlays(cOverlay).data(scanList);
  if isfield(overlays(cOverlay),'alphaOverlay') && ~isempty(overlays(cOverlay).alphaOverlay) && ~ismember(overlays(cOverlay).alphaOverlay,overlayNames)
    alphaOverlayName = overlays(cOverlay).alphaOverlay;
    cOverlay = cOverlay+1;
    overlayNames = [overlayNames {alphaOverlayName}];
    overlays = copyFields(viewGet(thisView,'overlay',viewGet(thisView,'overlaynum',alphaOverlayName)),overlays,cOverlay);
    overlays(cOverlay).data = overlays(cOverlay).data(scanList);
    disp(['(copyOverlay) Alpha overlay ''' alphaOverlayName ''' will also be copied']);
  end
end


if ~isempty(overlays)
  fprintf('(copyOverlay) %i overlays from %i scans have been copied\n',cOverlay,length(scanList));
else
  disp('(copyOverlay) No overlay have been copied');
end
return






