function roiBaseOverlayPlot(thisView,overlayNum,scanNum,x,y,z,roi)
%
%      
%        $Id$
%    authors: rs, jb 26/01/2010
%       goal: plot the distribution of current overlay values against current base anatomy values
%              uses box plot is the base anatomy only has two values (0 and 1)
%              use as an interrogate function in mrLoadRet
%

if isempty(thisView.ROIs)
   mrWarnDlg('(roiBaseOverlayPlot) Please load at least one ROI')
   return
end

% get underlying image intensity
baseNum = viewGet(thisView,'currentBase');
base = viewGet(thisView,'baseanatomy');
baseName = viewGet(thisView,'basename');
baseVoxelSize = viewGet(thisView,'baseVoxelSize');

%First, get a logical mask of the current overlay display, as well as the overlay data
nOverlays = length(overlayNum);
[mask,overlayData] = maskOverlay(thisView,overlayNum,scanNum);
mask = mask{1};
overlayData = overlayData{1};

%mask the overlay
overlayData(~mask) = NaN;

%get ROIs
roiList = viewGet(thisView,'visibleRois');
cRoi = 0;
for iRoi = roiList
  cRoi = cRoi+1;
  rois{cRoi} = viewGet(thisView,'roi',iRoi);
  base2roi(:,:,cRoi) = viewGet(thisView, 'base2roi',iRoi);
end
nRois = length(rois);

nTissueVoxels = zeros(nRois,nOverlays);
nVeinsVoxels = zeros(nRois,nOverlays);
nAllVoxels = zeros(nRois,nOverlays);
meanTissue = zeros(nRois,nOverlays);
meanVeins = zeros(nRois,nOverlays);


cOverlay=0;
for iOverlay = overlayNum
  overlayMax = 0;
  overlayMin = 0;
  cOverlay=cOverlay+1;
  %transform values in base space
  thisOverlayData{cOverlay} = getBaseSpaceOverlay(thisView, overlayData(:,:,:,cOverlay), scanNum, baseNum);
  dataSize = size(thisOverlayData{cOverlay});
  overlayName{cOverlay} = viewGet(thisView,'overlayName',iOverlay);
end

for iRoi = 1:nRois
  
  scancoords = xformROIcoords(rois{iRoi}.coords,inv(base2roi(:,:,iRoi)),rois{iRoi}.voxelSize,baseVoxelSize);
  % we need to remove any coordinate that might fall outside the base anatomy
  if ~isempty(scancoords)
    outside_voxels = find(scancoords(1,:)<1 | scancoords(1,:)>dataSize(1) |...
                      scancoords(2,:)<1 | scancoords(2,:)>dataSize(2) |...
                      scancoords(3,:)<1 | scancoords(3,:)>dataSize(3) );
    scancoords(:,outside_voxels) = [];
  end

  if ~isempty(scancoords)
    %RoiData
    scanRoiCoordsIndex = sub2ind(dataSize, scancoords(1,:)',  scancoords(2,:)',  scancoords(3,:)' );
    baseData = base.data(scanRoiCoordsIndex);
  end
  
  if ~isempty(scanRoiCoordsIndex)
    hFigure(iRoi) = figure('name',['ROI anatomy: ' rois{iRoi}.name]);
    
    for iOverlay = 1:cOverlay
      for jOverlay = iOverlay+1:cOverlay+1
        
        yData = thisOverlayData{iOverlay}(scanRoiCoordsIndex);
        yName = overlayName{iOverlay}(1:min(end,40));
        if jOverlay==cOverlay+1
          xData = baseData;
          xName=baseName;
        else
          xData = thisOverlayData{jOverlay}(scanRoiCoordsIndex);
          xName=overlayName{jOverlay};
        end
        xName = xName(1:min(end,50));
        nonNanValues = find(~isnan(yData)&~isnan(xData));
        yData = yData(nonNanValues);
        xData = xData(nonNanValues);

       hSubplot(iOverlay,jOverlay,iRoi) = subplot(nOverlays,nOverlays,(jOverlay-2)*nOverlays+iOverlay,'parent',hFigure(iRoi));
       hold on
        %see if base anatomy data are logical or continuous
        if length(unique(xData))<3 && jOverlay==cOverlay+1
           if length(unique(xData))==1
              vein_value = 1;
              tissue_value = 0;
           else
              vein_value = max(xData);
              tissue_value = min(xData);
           end
           scatter_plot = 0;

           %calculate number of voxel which are veins 
           vein_index = find(xData==vein_value);
           tissue_index = find(xData==tissue_value);
           labels = cell(size(yData));
           labels(tissue_index) = {'Tissue'};
           labels(vein_index) = {'Veins'};

           if ~isempty(vein_index) && ~isempty(tissue_index)
              boxplot(hSubplot(iOverlay,jOverlay,iRoi),yData,labels,'grouporder',{'Veins','Tissue'},'notch','on')
           end
           nTissueVoxels(iRoi,cOverlay) = length(tissue_index);
           nVeinsVoxels(iRoi,cOverlay) = length(vein_index);
           nAllVoxels(iRoi,cOverlay) = size(yData,1);

           title(hSubplot(iOverlay,jOverlay,iRoi),sprintf('%s (%d voxels, %.2f %% veins)', rois{iRoi}.name, nTissueVoxels(iRoi,cOverlay), nVeinsVoxels(iRoi,cOverlay)./nAllVoxels(iRoi,cOverlay)*100),'interpreter','none');

           %calculate mean and range of overlay for the veins and for the tissue
           %voxels
           fprintf(1,[rois{iRoi}.name '\n']);
           fprintf(1,['Veins: ' num2str(length(vein_index)) '/' num2str(nAllVoxels(iRoi,cOverlay)) ' voxels\t']);
           [string,meanVeins(iRoi,cOverlay)] = meanAndRange(yData(vein_index));
            strings{1} = ['Veins: ' string];
           fprintf(1,['Tissue: ' num2str(length(tissue_index)) '/' num2str(nAllVoxels(iRoi,cOverlay)) ' voxels\t']);
           [string,meanTissue(iRoi,cOverlay)] = meanAndRange(yData(tissue_index));
            strings{2} = ['Tissue: ' string];
            xlabel(hSubplot(iOverlay,jOverlay,iRoi),strings);

        else
           scatter_plot = 1;
           scatter(hSubplot(iOverlay,jOverlay,iRoi),xData,yData,2);
           %compute correlation
           [rho,p]=corrcoef([xData yData]);
           X=[ones(length(xData),1) xData];
           params = (X'*X)\X'*yData;
           if p(2)<.05
             xlim=get(hSubplot(iOverlay,jOverlay,iRoi),'xlim');
             ylim=get(hSubplot(iOverlay,jOverlay,iRoi),'ylim');
%              plot(hSubplot(iOverlay,jOverlay,iRoi),xData,X*params,'g');
             plot(hSubplot(iOverlay,jOverlay,iRoi),xlim,params(1)+params(2)*xlim,'r');
             set(hSubplot(iOverlay,jOverlay,iRoi),'ylim',ylim);
           end
           title(hSubplot(iOverlay,jOverlay,iRoi),sprintf('%d/%d voxels r=%.2f (p=%.3f) y=%.3f*x+%.3f',size(yData,1),length(scanRoiCoordsIndex),rho(2),p(2),params(2),params(1)),'interpreter','none');

           xlabel(hSubplot(iOverlay,jOverlay,iRoi),xName,'interpreter','none');
           ylabel(hSubplot(iOverlay,jOverlay,iRoi),yName,'interpreter','none');

%            base_max = max(base_max,max(baseData));
%            base_min = min(base_max,min(baseData));
        end
      end
%         overlayMax = max(overlayMax,max(roiOverlayData));
%         overlayMin = min(overlayMin,min(roiOverlayData));
    end
  end

%   for iRoi = 1:nRois
%      if ~isempty(nTissueVoxels(iRoi,cOverlay)) && ~isempty(nVeinsVoxels(iRoi,cOverlay))
%         subplot(hSubplot(cOverlay,iRoi))
%         scale = axis;
%         if scatter_plot
%            scale = [base_min base_max overlayMin overlayMax];
%         else
%            scale(3:4) = [overlayMin overlayMax];
%         end
%         axis(scale);
%      end
%   end
end
