% getNewSpaceOverlay.m
%
%      usage: newOverlayData = (overlayData, xform, newXCoords, newYCoords, newZCoords, interpMethod)
%
%         by: julien besle, loosely based on mrExport2SR by eli merriam
%       date: 11/05/2010
%    purpose: computes overlay values in a new space
%        $Id$	
%

function newOverlayData = getNewSpaceOverlay(overlayData, xform, newXCoords, newYCoords, newZCoords, interpMethod)

% interpMethod is used by calls to interp3 when
% extracting slices from the base and overlay arrays.
if ieNotDefined('interpMethod')
   interpMethod = mrGetPref('interpMethod');
end
if isempty(interpMethod)
    interpMethod = 'linear';
end

nOverlays = size(overlayData,4);
sliceDims(1) = size(newXCoords,1);
sliceDims(2) = size(newXCoords,2);
numPixels = sliceDims(1)*sliceDims(2);
newOverlayData = NaN([size(newXCoords) nOverlays]);

% round xform to the 7th decimal to avoid missing voxel due to round-off error
% when transformation should be a simple translation or the identity
epsilon = 1e-7;
xform = round(xform./epsilon).*epsilon;

hWaitbar = mrWaitBar(0,'Resampling overlay to new space');
% Compute new overlay data by base slice
nSlices = size(newXCoords,3);
scanDims = size(overlayData(:,:,:,1));
for iSlice = 1:nSlices
    
  baseCoordsHomogeneous = [reshape(newXCoords(:,:,iSlice),1,numPixels); reshape(newYCoords(:,:,iSlice),1,numPixels); reshape(newZCoords(:,:,iSlice),1,numPixels); ones(1,numPixels)];
  %baseCoords = reshape(baseCoordsHomogeneous(1:3,:)',[sliceDims 3]);

  % Transform coordinates
  overlayCoordsHomogeneous = xform * baseCoordsHomogeneous;
  %check if this slice has any data
  isInScan = overlayCoordsHomogeneous(1:3,:)>=0 & overlayCoordsHomogeneous(1:3,:)<=repmat(scanDims',1,size(overlayCoordsHomogeneous,2));
  
  if any(all(isInScan,1)) %if it does, interpolate the data form the overlay
    overlayCoords = reshape(overlayCoordsHomogeneous(1:3,:)',[sliceDims 3]);

    % Extract slice from current overlay.
    if ~isempty(overlayCoords) 
      for iOverlay=1:nOverlays
        newOverlayData(:,:,iSlice,iOverlay) = interp3(overlayData(:,:,:,iOverlay),overlayCoords(:,:,2),overlayCoords(:,:,1),overlayCoords(:,:,3), interpMethod,NaN);
        mrWaitBar((iSlice*nOverlays+iOverlay)/(nSlices*nOverlays),hWaitbar);
      end
    end
  end
end
mrCloseDlg(hWaitbar);
