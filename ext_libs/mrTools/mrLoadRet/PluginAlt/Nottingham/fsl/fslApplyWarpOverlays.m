function fslApplyWarpOverlays(thisView)
% applyWarpOverlays(thisView,thisView,overlayNum,scanNum,x,y,z)
%
%   applies non-linear FSL registration to overlay in current view and current analysis
%
% jb 23/07/2010
%
% $Id: applyWarpOverlays.m 2136 2011-05-31 10:23:12Z julien $ 

if strcmp(mrGetPref('fslPath'),'FSL not installed')
  mrWarnDlg('(fslApplyWarpOverlays) No path was provided for FSL. Please set MR preference ''fslPath'' by running mrSetPref(''fslPath'',''yourpath'')')
  return
end

keepAsking = 1;
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

[warpCoefFilename warpCoefPathname] = uigetfile('*.img;*.nii','FNIRT spline coefficient file');
if isnumeric(warpCoefFilename)
   return
end

%find temporary file extension based on FSL preference
switch(getenv('FSLOUTPUTTYPE'))
  case 'NIFTI'
    tempFilename='temp.nii';
  case 'NIFTI_PAIR'
    tempFilename='temp.img';
  case ''
    mrWarnDlg('(fslApplyWarpSurfOFF) Environment variable FSLOUTPUTTYPE is not set');
    return;
  otherwise
    mrWarnDlg(sprintf('(fslApplyWarpOverlays) Environment variable FSLOUTPUTTYPE is set to an unsupported value (%s)',getenv('FSLOUTPUTTYPE')));
    return;
end

warning('off','MATLAB:mat2cell:TrailingUnityVectorArgRemoved');
interpMethod = mrGetPref('interpMethod');
scanDims = viewGet(thisView,'scandims'); %we assume scans are all in the same space
scanHdr = viewGet(thisView,'niftiHdr');

overlayNames = viewGet(thisView,'overlayNames');
overlayNames = overlayNames(overlayList);
overlaysNumber = length(overlayList);
for iOverlay = 1:overlaysNumber
  overlays(iOverlay) = viewGet(thisView,'overlay',overlayList(iOverlay));
  if isfield(overlays(iOverlay),'alphaOverlay') && ~isempty(overlays(iOverlay).alphaOverlay) 
    alphaOverlayName = overlays(iOverlay).alphaOverlay;
    if ~ismember(overlays(iOverlay).alphaOverlay,overlayNames)
      overlaysNumber = overlaysNumber+1;
      overlayNames = [overlayNames {alphaOverlayName}];
      overlays(overlaysNumber) = viewGet(thisView,'overlay',viewGet(thisView,'overlaynum',alphaOverlayName));
      disp(['(fslApplyWarpOverlays) Alpha overlay ''' alphaOverlayName ''' will also be warped']);
    end
    overlays(iOverlay).alphaOverlay = ['warp(' alphaOverlayName ') - ' warpCoefFilename];
  end
end

for iOverlay = 1:overlaysNumber
   overlayIndices = false(1,length(scanList));
   for iScan = 1:length(scanList)
      overlayIndices(iScan) = ~isempty(overlays(iOverlay).data(scanList));
   end
   if any(overlayIndices)
      %transform the non empty overlay in a 4D matrix
      data = permute(reshape(cell2mat(overlays(iOverlay).data(scanList(overlayIndices))),[scanDims([1 2]) sum(overlayIndices) scanDims(3)]), [1 2 4 3]);
      warpedData = fslApplyWarp(data, [warpCoefPathname warpCoefFilename], tempFilename, scanHdr, interpMethod);
      overlays(iOverlay).data = cell(1,length(overlays(iOverlay).data));
      overlays(iOverlay).data(scanList(overlayIndices)) = mat2cell(warpedData,scanDims(1),scanDims(2),scanDims(3),ones(1,sum(overlayIndices)));
      overlays(iOverlay).name = ['warp(' overlays(iOverlay).name ') - ' warpCoefFilename];
      %this used to be helpful when it wasn't possible to readily identify which overlays are clipped
% % % %       %reset the clip values. This is due to the fact that warped data can go beyond the min and max of the original data (even when using nearest-neighboor, probably due to a bug of FSL
% % % %       overlays.clip = [min(warpedData(:)) max(warpedData(:))];

   end
end
thisView = viewSet(thisView,'newoverlay',overlays);
warning('on','MATLAB:mat2cell:TrailingUnityVectorArgRemoved');

refreshMLRDisplay(thisView.viewNum);




