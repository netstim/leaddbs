function computeImageSNR(thisView,overlayNum,scanNum,x,y,z,roi)
%
% computeImageSNR(thisView,thisView,overlayNum,scanNum,x,y,z)
%
%        $Id: findMaxRoi.m 2134 2011-05-31 08:49:17Z julien $
% jb 01/06/2011
% Computes image SNR as the ratio of the signal mean in a gray matter ROI and
%     the signal std deviation across voxels in a white matter ROI 
%     iSNR is then averaged across all samples of the scan(s).

scanList = [1:viewGet(thisView,'numscans')];

whiteRoiNum = [];
greyRoiNum = [];
keepAsking=true;

while keepAsking
  whiteRoiNum = selectInList(thisView,'rois','Select White Matter ROI',whiteRoiNum);
  while length(whiteRoiNum)>1
    mrWarnDlg('(computeImageSNR) Please select only one ROI');
    whiteRoiNum = selectInList(thisView,'rois','Select White Matter ROI',whiteRoiNum);
  end
  if isempty(whiteRoiNum)
    return;
  end
  
  while keepAsking
    greyRoiNum = selectInList(thisView,'rois','Select Grey Matter ROI',greyRoiNum);
    while length(greyRoiNum)>1
      mrWarnDlg('(computeImageSNR) Please select only one ROI');
      greyRoiNum = selectInList(thisView,'rois','Select Grey Matter ROI',greyRoiNum);
    end
    if isempty(greyRoiNum) && size(greyRoiNum,2)==1
      return;
    elseif size(greyRoiNum,2)==0
      break;
    end
    
    while keepAsking
      if length(scanList)>1
        scanList = selectInList(thisView,'scans',[],scanList);
        if isempty(greyRoiNum) && size(greyRoiNum,2)==1
          return;
        elseif isempty(scanList)
          break;
        end
      end
      keepAsking=false;
    end
  end
end

[d,whichRoi,roiCoords] = loadScanRois(thisView,scanList,[whiteRoiNum greyRoiNum]);
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

  % get the indices of the voxels of this ROI that are not masked
  unmaskedWhiteVoxels = mask(sub2ind(size(mask),roiCoords{1}(1,:)',roiCoords{1}(2,:)',roiCoords{1}(3,:)'));
  whiteData = permute(d.data(whichRoi{1}(unmaskedWhiteVoxels),:,:,:),[1 4 2 3]);
  unmaskedGreyVoxels = mask(sub2ind(size(mask),roiCoords{2}(1,:)',roiCoords{2}(2,:)',roiCoords{2}(3,:)'));
  greyData = permute(d.data(whichRoi{2}(unmaskedGreyVoxels),:,:,:),[1 4 2 3]);
else
  whiteData = permute(d.data(whichRoi{1},:,:,:),[1 4 2 3]);
  greyData = permute(d.data(whichRoi{2},:,:,:),[1 4 2 3]);
end
  
whiteRoiNumberValues = nnz(all(~isnan(whiteData),2));
 whiteRoiNumberVoxels = nnz(whichRoi{1});

greyRoiNumberValues = nnz(all(~isnan(greyData),2));
greyRoiNumberVoxels = nnz(whichRoi{2});

% meanWhite = nanmean(whiteData(:));
% stdWhite = mean(nanstd(whiteData,0));
% meanGrey =  mean(greyData(:));
% stdGrey =  mean(nanstd(greyData,0));
% imageSNR = mean(nanmean(greyData)./nanstd(whiteData,0));
% imageSNRstd = std(nanmean(greyData)./nanstd(whiteData,0));
  
fprintf(1,['\nROI White Matter: ' viewGet(thisView,'roiName',whiteRoiNum) '(' num2str(whiteRoiNumberValues) '/' num2str(whiteRoiNumberVoxels) ' scan voxels)\n']);
fprintf(1,['ROI Grey Matter: ' viewGet(thisView,'roiName',greyRoiNum) '(' num2str(greyRoiNumberValues) '/' num2str(greyRoiNumberVoxels) ' scan voxels)\n']);
% fprintf(1,['\t\t White matter signal average :' num2str(meanWhite) ' (average std-dev across samples = ' num2str(stdWhite) ')\n']);
% fprintf(1,['\t\t Grey matter signal average :' num2str(meanGrey) ' (average std-dev across samples = ' num2str(stdGrey) ')\n']);
% fprintf(1,['\t\tAverage image SNR :' num2str(imageSNR) ' (std dev across samples = ' num2str(imageSNRstd) ')\n']);
  

fprintf(1,'\tROI average (+/- stddev), averaged across frames \t image SNR (average +/- stddev across frames)\n');
fprintf(1,'\t\t\tWhite matter\t\tGrey matter\n');

cScan=0;
for iScan=scanList
  cScan=cScan+1;
  thisWhiteData = whiteData(:,runs(cScan,1):runs(cScan,2));
  thisGreyData = greyData(:,runs(cScan,1):runs(cScan,2));
  meanWhite(cScan) = nanmean(thisWhiteData(:));
  stdWhite(cScan) = mean(nanstd(thisWhiteData,0));
  meanGrey(cScan) =  mean(thisGreyData(:));
  stdGrey(cScan) =  mean(nanstd(thisGreyData,0));
  imageSNR(cScan) = mean(nanmean(thisGreyData)./nanstd(thisWhiteData,0));
  imageSNRstd(cScan) = std(nanmean(thisGreyData)./nanstd(thisWhiteData,0));
  fprintf(['\t  Scan ' num2str(iScan) ':\t' displayValues(meanWhite(cScan),stdWhite(cScan),meanGrey(cScan),stdGrey(cScan),imageSNR(cScan),imageSNRstd(cScan)) '\n']);
end

if length(scanList)>1
  fprintf('\tAverage across scans (weighted by scan length):\n');
  meanWhite = sum(nFrames/runs(end).*meanWhite);
  stdWhite = sum(nFrames/runs(end).*stdWhite);
  meanGrey = sum(nFrames/runs(end).*meanGrey);
  stdGrey = sum(nFrames/runs(end).*stdGrey);
  imageSNR = sum(nFrames/runs(end).*imageSNR);
  imageSNRstd = sum(nFrames/runs(end).*imageSNRstd);
  fprintf('\t\t\t%s\n\n',displayValues(meanWhite,stdWhite,meanGrey,stdGrey,imageSNR,imageSNRstd));
end



function string = displayValues(average1,stddev1,average2,stddev2,tSNR,tSNRstd)

string = sprintf('%.0f (+/- %.0f)\t%.0f (+/- %.0f)\t %.1f (+/- %.2f)',average1,stddev1,average2,stddev2,tSNR,tSNRstd);

