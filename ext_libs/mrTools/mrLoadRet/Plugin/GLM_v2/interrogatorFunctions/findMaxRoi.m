function findMaxRoi(thisView,overlayNum,scanNum,x,y,z,roi)
%
% findMaxContiguousVoxels(thisView,thisView,overlayNum,scanNum,x,y,z)
%
%        $Id$
% jb 04/12/2009
% finds the max and min value, as well as coordinates of all non-zero voxels in  ROI
% (non-zero meaning between any boundaries that have been chosen from any
% overlay)  
%
if isempty(roi)
   mrWarnDlg('(findMaxRoi) Please click in an ROI')
   return
end


baseNum = viewGet(thisView,'currentBase');

if isempty(overlayNum)
   mrWarnDlg('(findMaxRoi) Please load an overlay')
   return
end

%First, get a logical mask of the current overlay display, as well as the overlay data
[mask,overlayData] = maskOverlay(thisView,overlayNum,scanNum);


mask = mask{1};
overlayData = overlayData{1};
%mask the overlay
overlayData(~mask) = NaN;

base2scan = viewGet(thisView,'base2scan',scanNum,[],baseNum);
base2tal =  viewGet(thisView,'base2tal',baseNum);
if ~isempty(base2tal)
   scan2tal = base2tal/base2scan;
end
%scan2mag = viewGet(thisView,'scanXform',scanNum);
scanvoxelsize = viewGet(thisView,'scanvoxelsize');


cOverlay=0;
for iOverlay = overlayNum;
  cOverlay= cOverlay+1;
  overlayName = viewGet(thisView,'overlayName',iOverlay);
  thisOverlayData = overlayData(:,:,:,cOverlay);

  disp(['Overlay: ' overlayName])

  for iRoi = 1:length(roi)

     scanRoiCoords = getROICoordinates(thisView,viewGet(thisView,'roinum',roi{iRoi}.name));
     if ~isempty(scanRoiCoords)   
        scanRoiCoordsIndex = sub2ind(size(thisOverlayData), scanRoiCoords(1,:)',  scanRoiCoords(2,:)',  scanRoiCoords(3,:)' );
        roiOverlayData = thisOverlayData(scanRoiCoordsIndex);
        totalNumberVoxels = numel(roiOverlayData);

        [max_value max_index] = max(roiOverlayData);
        [min_value min_index] = min(roiOverlayData);
        max_coordinates = scanRoiCoords(:,max_index);
        min_coordinates = scanRoiCoords(:,min_index);

        %find the center of mass
        scanRoiCoordsIndex = scanRoiCoordsIndex(~isnan(roiOverlayData));
        roiOverlayData = roiOverlayData(~isnan(roiOverlayData));
        [coordsX,coordsY,coordsZ] = ind2sub(size(thisOverlayData),scanRoiCoordsIndex);
        COM_coordinates = round(sum(repmat(roiOverlayData,1,3)/sum(roiOverlayData) .* [coordsX,coordsY,coordsZ],1))';

        max_base_coordinates = (base2scan\[max_coordinates;1]);
        min_base_coordinates = (base2scan\[min_coordinates;1]);
        COM_base_coordinates = (base2scan\[COM_coordinates;1]);
        

        fprintf(1,['\tROI ' roi{iRoi}.name '(' num2str(numel(roiOverlayData)) '/' num2str(totalNumberVoxels) ' scan voxels)\n']);
        fprintf(1,['\t\tmax value :' num2str(max_value) '\n']);
        fprintf(1,['\t\tmax scan coordinates :' num2str(max_coordinates(1:3)') '\n']);
        fprintf(1,['\t\tmax base coordinates :' num2str(round(max_base_coordinates(1:3)')) '\n']);
        if ~isempty(base2tal)
           max_tal_coords = round(scan2tal*max_base_coordinates);
           fprintf(1,['\t\tmax Talairach coordinates: ' num2str(max_tal_coords(1:3)') '\n']);
        end

        fprintf(1,['\t\tmin value :' num2str(min_value) '\n']);
        fprintf(1,['\t\tmin scan coordinates :' num2str(min_coordinates(1:3)') '\n']);
        fprintf(1,['\t\tmin base coordinates :' num2str(round(min_base_coordinates(1:3)')) '\n']);
        if ~isempty(base2tal)
           min_tal_coords = round(scan2tal*min_base_coordinates);
           fprintf(1,['\t\tmin Talairach coordinates: ' num2str(min_tal_coords(1:3)') '\n']);
        end
        fprintf(1,'\n');

        fprintf(1,['\t\tcenter-of-mass scan coordinates :' num2str(COM_coordinates(1:3)') '\n']);
        fprintf(1,['\t\tcenter-of-mass base coordinates :' num2str(round(COM_base_coordinates(1:3)')) '\n']);
        if ~isempty(base2tal)
           COM_tal_coords = round(scan2tal*COM_base_coordinates);
           fprintf(1,['\t\tcenter-of-mass Talairach coordinates: ' num2str(COM_tal_coords(1:3)') '\n']);
        end
        fprintf(1,'\n');

        fprintf(1,['\t\taverage across ROI :' num2str(nanmean(roiOverlayData)) '\n']);
        fprintf(1,['\t\tstandard deviation :' num2str(nanstd(roiOverlayData)) '\n']);
        fprintf(1,['\t\trms across ROI :' num2str(sqrt(nanmean(roiOverlayData.^2))) '\n']);
        fprintf(1,'\n');

     else
        fprintf(1,[roi{iRoi}.name ' is empty\n\n']);
     end
  end
  
end






