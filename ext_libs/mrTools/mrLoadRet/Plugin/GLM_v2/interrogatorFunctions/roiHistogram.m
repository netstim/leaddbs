function roiHistogram(thisView,overlayNum,scanNum,x,y,z,roi)
%
% roiHistogram(thisView,thisView,overlayNum,scanNum,x,y,z)
%
%        $Id$
% jb 10/12/2009
% plots the histogram of all non-NaN overlay values across voxels of all
% visible ROIs
% (non-NaN meaning between any boundaries that have been chosen from any
% overlay) 
%

if length(overlayNum)>1
  mrWarnDlg('(roiHistogram) Not implemented for several overlays');
  return;
end
  
  
% Error if no current ROI
if isempty(roi)
	disp('Voxel outside an ROI, using all visible ROIs');
   n_rois =  viewGet(thisView,'numberofROIs');
   if n_rois
      roi{1} = viewGet(thisView,'roi',1);
      if n_rois >1
         roi{1}.name = 'All ROIs';
         for i_roi = 2:n_rois
            temp_roi = viewGet(thisView,'roi',i_roi);
            roi{1}.coords = union(roi{1}.coords',temp_roi.coords','rows')';
            %here should check that ROIs are compatible
         end
      end
   else
      mrWarnDlg('(roiHistogram) No ROI currently loaded.');
   end
end

% set roi coords
for roinum = 1:length(roi)
  % get scanNum coordinates
  roi{roinum}.scanCoords = getROICoordinates(thisView,roi{roinum},scanNum);
end

if isempty(overlayNum)
   mrWarnDlg('(roiHistogram) Please load an overlay')
   return
end
[mask,overlayData] = maskOverlay(thisView,overlayNum,scanNum);
mask=mask{1};
overlayData = overlayData{1};

overlayNames = viewGet(thisView,'overlayNames');
[mask_coordinates(:,1) mask_coordinates(:,2) mask_coordinates(:,3)] = ind2sub(size(mask),find(mask));


% select the window to plot into
fignum = selectGraphWin;

% turn off menu/title etc.
set(fignum,'NumberTitle','off');
set(fignum,'Name','roiHistogram');

for i_roi = 1:length(roi)
  subplot(1,length(roi),i_roi);
  coords = intersect(roi{i_roi}.scanCoords',mask_coordinates,'rows');
  overlay_values = overlayData(sub2ind(size(overlayData),coords(:,1), coords(:,2), coords(:,3)));
  
  hist(overlay_values,50);
  if strcmp(overlayNames{overlayNum},'ph')
     %rose(overlay_values,50);
     scale = axis;
     scale(1:2) = [0 2*pi];
     axis(scale);
     set(gca,'XTick',0:pi/2:2*pi)
     set(gca,'XTickLabel',{'0','pi/2','pi','3pi/2','2pi'})
  end
  xlabel(overlayNames{overlayNum});
  ylabel('Number of voxels');

  title([roi{i_roi}.name ' (' num2str(size(coords,1)) '/' num2str(size(roi{i_roi}.scanCoords,2)) ' voxels)']);
  
  %print values
  disp(overlayNames{overlayNum});
  [binCount,binCenter] = hist(overlay_values,50);
  binCenter = binCenter(binCount>0);
  binCount = binCount(binCount>0);
  disp(binCenter');
  disp(binCount');
end





