function findMaxRoi(thisView,overlayNum,scanNum,x,y,z,roi)
%
% getMaxRoiHrf(thisView,thisView,overlayNum,scanNum,x,y,z)
%
%        $Id: findMaxRoi.m 2058 2011-02-24 16:44:39Z julien $
% jb 04/12/2009
% finds the max and min value, as well as coordinates of all non-zero voxels in  ROI
% (non-zero meaning between any boundaries that have been chosen from any
% overlay)  
%
if isempty(roi)
   mrWarnDlg('(findMaxRoi) Please click in an ROI')
   return
end

analysisType = viewGet(thisView,'analysisType');
glmParams = convertOldGlmParams(viewGet(thisView,'analysisParams'));

if ~ismember(analysisType,{'glmAnalStats','glmAnal','glmcAnal','erAnal','deconvAnal'})
   disp(['(glmPlot) Wrong type of analysis (' analysisType ')']);
   return;
end

dt=.05;
d = viewGet(thisView,'d');
[~,hrf] = feval(glmParams.hrfModel,glmParams.hrfParams,dt,0,1);
scanDims = viewGet(thisView,'scandims');

betas = permute(d.ehdr,[4 5 1 2 3]);
for iRoi=0:length(roi)
  if iRoi
    fprintf('\nRoi %s:\n',roi{iRoi}.name);
    roi{iRoi}.scanCoords = getROICoordinates(thisView,roi{iRoi},scanNum);
    indices = roi{iRoi}.scanCoords(1:3,:)';
  else
    fprintf('\nVoxel %i %i %i:\n',x,y,z);
    indices=[x y z];
  end
  indices = sub2ind(scanDims,indices(:,1),indices(:,2),indices(:,3));
  
  meanHrf = hrf*nanmean(betas(:,:,indices),3)';
  [maxRoiHrf,maxSample] = max(meanHrf);
  [minRoiHrf,minSample] = min(meanHrf);
  fprintf('\t\t\tMax Hrf:');
  fprintf('%.4f ',maxRoiHrf);
  fprintf('\n\tLatency Max Hrf:');
  fprintf('%.2f   ',dt*maxSample);
  fprintf('\n');
  fprintf('\t\t\tMin Hrf:');
  fprintf('%.4f ',minRoiHrf);
  fprintf('\n\tLatency Min Hrf:');
  fprintf('%.2f   ',dt*minSample);
  fprintf('\n');
  
end
figure;plot(dt*[1:size(meanHrf,1)],meanHrf);


