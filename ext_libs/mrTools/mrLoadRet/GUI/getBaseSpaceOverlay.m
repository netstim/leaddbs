% getBaseSpaceOverlay.m
%
%      usage: [newOverlayData, baseVoxelSize, baseCoords] = getBaseSpaceOverlay(thisView,overlayData,scanNum,baseNum,interpMethod, depthBins, rotateAngle)
%         by: julien besle, loosely based on mrExport2SR by eli merriam
%       date: 25/01/2010
%    purpose: computes overlay values in a base anatomy space
%        $Id$ 
%     output: newOverlayData: interpolated values
%             baseVoxelSize: size of voxels in base anatomy (estimated in the case of flat maps)
%             baseCoords: base coordinates arranged in the same volume as the returned data
%                         (useful to remap interpolated values into 3D space for flat maps)
%

function [newOverlayData, baseVoxelSize, baseCoords] = getBaseSpaceOverlay(thisView,overlayData,scanNum,baseNum,interpMethod, depthBins, rotateAngle)

if ieNotDefined('scanNum')
    scanNum = viewGet(thisView,'curscan');
end
if ieNotDefined('baseNum')
    baseNum = viewGet(thisView,'curbase');
end
% interpMethod and interpExtrapVal are used by calls to interp3 when
% extracting slices from the base and overlay arrays.
if ieNotDefined('interpMethod')
   interpMethod = mrGetPref('interpMethod');
end
if ieNotDefined('interpMethod')
    interpMethod = 'linear';
end

basedims = viewGet(thisView, 'basedims', baseNum);
base2scan = viewGet(thisView,'base2scan',scanNum,[],baseNum);
baseType = viewGet(thisView,'baseType',baseNum);

baseVoxelSize = viewGet(thisView,'basevoxelsize',baseNum);

if all(all(abs(base2scan - eye(4))<1e-6)) %check if we're in the scan space
   disp('(getBaseSpaceOverlay) Overlay already in the base space, skipping resampling...');
   newOverlayData = overlayData;
   return;
end


switch(baseType)
   case 0   %the base is a volume
   % Generate coordinates with meshgrid
   [Ycoords,Xcoords,Zcoords] = meshgrid(1:basedims(2),1:basedims(1),1:basedims(3));
   
   case 1   %the base is a flat map
      sliceIndex = viewGet(thisView,'baseSliceIndex',baseNum);
      if ieNotDefined('depthBins')
         depthBins = viewGet(thisView,'corticalDepthBins');
      end
      if ieNotDefined('rotateAngle')
         rotateAngle = viewGet(thisView,'rotate');
      end
      basedims(3) = depthBins;
      Xcoords = zeros(basedims(1),basedims(2),basedims(3));
      Ycoords = Xcoords; Zcoords=Xcoords;
      baseCoordMap = viewGet(thisView,'baseCoordMap',baseNum,0:1/(depthBins-1):1);
      if ~isempty(baseCoordMap)
        % only use the baseCoordMap for when the slice
        % in the third dimension (no other view of a 
        % flat map is really valid).
        if sliceIndex == 3
          Xcoords0 = permute(baseCoordMap.coords(:,:,1,1,:),[1 2 5 3 4]);
          Ycoords0 = permute(baseCoordMap.coords(:,:,1,2,:),[1 2 5 3 4]);
          Zcoords0 = permute(baseCoordMap.coords(:,:,1,3,:),[1 2 5 3 4]);
        else
          oneTimeWarning('badSliceIndex',sprintf('(getBaseSpaceOverlay) Trying to display a flat/surface with the sliceIndex set to %i instead of 3. This is probably because there is something wrong with the Nifti Qforms at your site -- specifically you should check in mrAlign whether your volume displays correctly for when you have click the saggital, coronal and axial buttons. If not, you will need to swap dimensions until they do and then make sure all of your qforms have their dimensions in the same order. Your overlays will not display correctly on this volume.',sliceIndex));
          %attempt to do it anyway
          try
            Xcoords0 = permute(baseCoordMap.coords(:,:,1,1,:),[1 2 5 3 4]);
            Ycoords0 = permute(baseCoordMap.coords(:,:,1,2,:),[1 2 5 3 4]);
            Zcoords0 = permute(baseCoordMap.coords(:,:,1,3,:),[1 2 5 3 4]);
          catch errorId
            mrErrorDlg(errorId.message)
          end
        end
        %rotate coordinates
        if  rotateAngle
          for iDepth = 1:depthBins
             Xcoords(:,:,iDepth) = imrotate(Xcoords0(:,:,iDepth),rotateAngle,'bilinear','crop');
             Ycoords(:,:,iDepth) = imrotate(Ycoords0(:,:,iDepth),rotateAngle,'bilinear','crop');
             Zcoords(:,:,iDepth) = imrotate(Zcoords0(:,:,iDepth),rotateAngle,'bilinear','crop');
          end
        else
          Xcoords=Xcoords0;
          Ycoords=Ycoords0;
          Zcoords=Zcoords0;
        end
      end
      
      %just as an indication, the voxel size is the mean distance between voxels consecutive in the 3 directions
      %mask coordinates with non-rotated coordinates (to avoid edges introduced by imrotate)
      Xcoords0Mask = Xcoords0==0;
      Xcoords0Mask = convn(Xcoords0Mask,ones(5,5,5),'same'); %expand the mask a bit to make sure we don't include any edge voxels
      XcoordsNaN = Xcoords; 
      XcoordsNaN(Xcoords0Mask>0)=NaN;
      YcoordsNaN = Ycoords;
      YcoordsNaN(Xcoords0Mask>0)=NaN;
      ZcoordsNaN = Zcoords;
      ZcoordsNaN(Xcoords0Mask>0)=NaN;
      baseVoxelSize(1) = baseVoxelSize(1)*nanmean(nanmean(nanmean(sqrt(diff(XcoordsNaN,1,1).^2 + diff(YcoordsNaN,1,1).^2 + diff(ZcoordsNaN,1,1).^2))));
      baseVoxelSize(2) = baseVoxelSize(2)*nanmean(nanmean(nanmean(sqrt(diff(XcoordsNaN,1,2).^2 + diff(YcoordsNaN,1,2).^2 + diff(ZcoordsNaN,1,2).^2))));
      baseVoxelSize(3) = baseVoxelSize(3)*nanmean(nanmean(nanmean(sqrt(diff(XcoordsNaN,1,3).^2 + diff(YcoordsNaN,1,3).^2 + diff(ZcoordsNaN,1,3).^2))));

   otherwise
      mrWarnDlg('(getBaseSpaceOverlay) This function is not implemented for surfaces')
      
end

newOverlayData = getNewSpaceOverlay(overlayData, base2scan, Xcoords, Ycoords, Zcoords, interpMethod);

if nargout==3
   baseCoords = cat(4,Xcoords, Ycoords);
   baseCoords = cat(4,baseCoords, Zcoords);
end
   
   


