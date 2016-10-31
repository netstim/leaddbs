function computeTemporalSNR(thisView,overlayNum,scanNum,x,y,z,roi)
%
% computeTemporalSNR(thisView,thisView,overlayNum,scanNum,x,y,z)
%
%        $Id: findMaxRoi.m 2134 2011-05-31 08:49:17Z julien $
% jb 01/06/2011
% Computes temporal SNR as the ratio of the signal mean across time samples
%     to the signal std deviation across time samples of the selected scan(s)
%     tSNR is then averaged across all voxels of each ROI
%

scanList = 1:viewGet(thisView,'numScans');

% version where user has to click in the ROI
% % if isempty(roi)
% %    mrWarnDlg('(computeTemporalSNR) Please click in an ROI')
% %    return
% % end
% % for iRoi= 1:length(roi)
% %   roiList(iRoi) = viewGet(thisView,'roiNum',roi{iRoi}.name);
% % end
% % 
% % if length(scanList)>1
% %   scanList = selectInList(thisView,'scans',[],scanList);
% %   if isempty(scanList)
% %     return;
% %   end
% % end

% version where user chooses ROIs in a list
keepAsking=true;
roiList = [];
while keepAsking
  roiList = selectInList(thisView,'rois','Select ROI(s)',roiList);
  if isempty(roiList)
    return;
  end
  while keepAsking
    if length(scanList)>1
      scanList = selectInList(thisView,'scans',[],scanList);
      if isempty(scanList)
        break;
      end
    end
    keepAsking=false;
  end
end

[d,whichRoi,roiCoords] = loadScanRois(thisView,scanList,roiList);
cScan=0;
endFrame=0;
for iScan=scanList
  cScan = cScan+1;
  nFrames(cScan) = viewGet(thisView,'nframes',iScan);
  runs(cScan,:) = endFrame+ [1 nFrames(cScan)];
  endFrame = endFrame+nFrames(cScan);
end


if ~isempty(overlayNum)
  %get the current mask
  mask = maskOverlay(thisView,overlayNum,scanNum);
  mask = mask{1};
  if all(~mask(:)) %if everything is masked, then nothing is masked
    mask(:) = 1;
  end
end

for iRoi=1:length(roiList)
  roiNumberVoxels = nnz(whichRoi{iRoi});
  if ~isempty(overlayNum)
    % get the indices of the voxels of this ROI that are not masked
    unmaskedVoxels = mask(sub2ind(size(mask),roiCoords{iRoi}(1,:)',roiCoords{iRoi}(2,:)',roiCoords{iRoi}(3,:)'));
    data = permute(d.data(whichRoi{iRoi}(unmaskedVoxels),:,:,:),[1 4 2 3]);
  else
    data = permute(d.data(whichRoi{iRoi},:,:,:),[1 4 2 3]);
  end
  %data = permute(d.data(whichRoi{iRoi},:,:,:),[1 4 2 3]);
  roiNumberValues = nnz(all(~isnan(data),2));
  
  fprintf(1,['\nROI ' viewGet(thisView,'roiName',roiList(iRoi)) '(' num2str(roiNumberValues) '/' num2str(roiNumberVoxels) ' scan voxels)\n']);
  fprintf(1,'\tTemporal average (+/- stddev), averaged across voxels \t Temporal SNR (average +/- stddev across voxels)\n');

  cScan=0;
  for iScan=scanList
    cScan=cScan+1;
    thisData = data(:,runs(cScan,1):runs(cScan,2));
    meanScanSignal(cScan) = nanmean(thisData(:));
    stdScanSignal(cScan) = mean(nanstd(thisData,0,2));
    tScanSNR(cScan) = mean(nanmean(thisData,2)./nanstd(thisData,0,2));
    tScanSNRstd(cScan) = std(nanmean(thisData,2)./nanstd(thisData,0,2));
    fprintf(['\t  Scan ' num2str(iScan) ':\t' displayValues(meanScanSignal(cScan),stdScanSignal(cScan),tScanSNR(cScan),tScanSNRstd(cScan)) '\n']);
  end
 
  if length(scanList)>1
    fprintf('\tAverage across scans (weighted by scan length)\n');
    meanSignal = sum(nFrames/runs(end).*meanScanSignal);
    stdSignal = sum(nFrames/runs(end).*stdScanSignal);
    tSNR = sum(nFrames/runs(end).*tScanSNR);
    tSNRstd = sum(nFrames/runs(end).*tScanSNRstd);
    fprintf('\t\t\t%s\n\n',displayValues(meanSignal,stdSignal,tSNR,tSNRstd));
  end
end


function string = displayValues(average,stddev,tSNR,tSNRstd)

string = sprintf('%.0f (+/- %.0f)\t\t\t%.1f (+/- %.1f)',average,stddev,tSNR,tSNRstd);






