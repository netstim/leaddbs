function eventRelatedContrast(view,overlayNum,scan,x,y,s,roi)
%
%      usage: 
%         by: Luke Hallum
%       date: 06/06/2008
%    purpose: Interrogator that plots derived HRFs, their 
%             contrast, and assess the statistical significance
%             of that contrast via a resampling method.

mrGlobals

CONDCONTRAST = [1 2];
NUMRESAMPLES = 1000;
EPOCH = 20;
ALPHA = 0.05;
[B,A] = butter(5,[0.05 0.9]);
curGroup = viewGet(view,'curgroup');
curScan = viewGet(view,'curscan');
tr_s = MLR.groups(curGroup).scanParams(curScan).framePeriod;
vectorTime = tr_s * ((1:EPOCH) - 0.5); % for plotting

% check arguments
if ~any(nargin == [1:7])
  help eventRelatedContrast
  return;
end

% select the window to plot into
fignum = selectGraphWin;

% turn off menu/title etc.
set(fignum,'NumberTitle','off');
set(fignum,'Name','eventRelatedContrast');

% set roi coords
for roinum = 1:length(roi)
  % get scan coordinates
  roi{end}.scanCoords = getROICoordinates(view,roi{roinum},scan);
end

%%%%%%%%%%
% Get time series data for this voxel, and all voxels w/i current 
% ROI.
%%%%%%%%%%
datafile = viewGet(view,'tseriespathstr',curScan,curGroup);
coords = roi{1}.scanCoords;
handleWaitBar = mrWaitBar(0,'Loading data...');
for ii = 1:size(coords,2)

  [dd,hh] = mlrImageReadNifti(datafile,{[coords(1,ii)],[coords(2,ii)],[coords(3,ii)],[]},'double',1);
  voxelTS(:,ii) = squeeze(dd);

  mrWaitBar(ii/size(coords,2),handleWaitBar);
end
close(handleWaitBar);
[dd,hh] = mlrImageReadNifti(datafile,{[x],[y],[s],[]},'double',1);
voxelTSthis = squeeze(dd);
%
%%%%%%%%%%

%%%%%%%%%%
% Build design matrix from the 'stim' files, then build 
% convolution matrix.
%%%%%%%%%%
junkedFrames = viewGet(view,'totaljunkedframes',curScan,curGroup);
stimfileData = viewGet(view,'stimfile',curScan,curGroup);
numConditions = size(stimfileData{1}.mylog.stimtimes_s,2);
lenScanLessJunk = size(voxelTS,1)/size(stimfileData,2);
dsgM = zeros(size(voxelTS,1),numConditions);
for ii = 1:size(stimfileData,2)

  for kk = 1:size(dsgM,2)

    aa = stimfileData{ii}.mylog(1).stimtimes_s{kk};
    jj = find(aa > junkedFrames(ii) * tr_s);
    dsgM(round(aa(jj)/tr_s - junkedFrames(ii) + (ii-1)*lenScanLessJunk),kk) = 1;
  end
end
dsgMconvmtx = erConGenConvMtx(dsgM,EPOCH);
%
%%%%%%%%%%

%%%%%%%%%%
% Solve for HRFs; plot those HRFs and their contrasts.
%%%%%%%%%%
voxelTSthis = filtfilt(B,A,voxelTSthis);
voxelTS = filtfilt(B,A,voxelTS);
thetathis = dsgMconvmtx \ voxelTSthis;
thetathis = reshape(thetathis,EPOCH,numConditions);
theta = mean(dsgMconvmtx \ voxelTS,2);
theta = reshape(theta,EPOCH,numConditions);
subplot(3,2,1); plot(vectorTime,100*thetathis(:,CONDCONTRAST),'LineWidth',2.0); title(sprintf('Derived HRFs voxel %d, %d, %d',x,y,s)); ylabel(sprintf('fMRI signal [per cent signal change]')); drawnow
subplot(3,2,2); plot(vectorTime,100*theta(:,CONDCONTRAST),'LineWidth',2.0); title(sprintf('Derived HRFs ROI "%s"',roi{1}.name)); drawnow
subplot(3,2,3); plot(vectorTime,100*(thetathis(:,CONDCONTRAST(2)) - thetathis(:,CONDCONTRAST(1))),'-k','LineWidth',2.0); hold on; title(sprintf('Contrast of derived HRFs for voxel')); drawnow
subplot(3,2,4); plot(vectorTime,100*(theta(:,CONDCONTRAST(2)) - theta(:,CONDCONTRAST(1))),'-k','LineWidth',2.0); hold on; title(sprintf('Contrast of derived HRFs for ROI')); drawnow
%
%%%%%%%%%%

%%%%%%%%%%
% Here, we repeat the above solution many times, each time using
% a different design matrix, generated at random. That is, we're
% resampling the design matrix (and effectively resampling the
% contrasts of interest).
%%%%%%%%%%
handleWaitBar = mrWaitBar(0,'Resampling contrasts...');
for kk = 1:NUMRESAMPLES

  dsgM = erConRandDsg1(dsgM);
  dsgMconvmtx = erConGenConvMtx(dsgM,EPOCH);

  resampleThetaThis(:,kk) = dsgMconvmtx \ voxelTSthis;
  aa = reshape(resampleThetaThis(:,kk),EPOCH,numConditions);
  contrastThis(:,kk) = aa(:,CONDCONTRAST(2)) - aa(:,CONDCONTRAST(1));

  resampleTheta(:,kk) = mean(dsgMconvmtx \ voxelTS,2);
  aa = reshape(resampleTheta(:,kk),EPOCH,numConditions);
  contrast(:,kk) = aa(:,CONDCONTRAST(2)) - aa(:,CONDCONTRAST(1));

  mrWaitBar(kk/NUMRESAMPLES,handleWaitBar);
end
close(handleWaitBar);
aa = reshape(mean(resampleThetaThis,2),EPOCH,numConditions);
subplot(3,2,5); plot(vectorTime,100*contrastThis,'-k','LineWidth',0.5); hold on; title(sprintf('Resampled contrasts for voxel')); xlabel(sprintf('Time from probe onset [s]'))
plot(vectorTime,100*mean(contrastThis,2),'-r','LineWidth',2.0);
aa = reshape(mean(resampleTheta,2),EPOCH,numConditions);
subplot(3,2,6); plot(vectorTime,100*contrast,'-k','LineWidth',0.5); hold on; title(sprintf('Resampled contrasts for ROI'));
plot(vectorTime,100*mean(contrast,2),'-r','LineWidth',2.0);
subplot(3,2,4); plot(vectorTime,100*mean(contrast,2),'-r','LineWidth',2.0);
subplot(3,2,3); plot(vectorTime,100*mean(contrastThis,2),'-r','LineWidth',2.0);
%
%%%%%%%%%%

%%%%%%%%%%
% For each estimated parameter (that is, the values of HRFs at
% each point in time), locate the tails that contain only
% ALPHA*100 per cent of the mass of the resampled distribution.
%%%%%%%%%%
interval95mass = zeros(size(contrast,1),2);
interval95massThis = zeros(size(contrastThis,1),2);
for ii = 1:size(contrast,1)

  [count,binCentre] = hist(contrast(ii,:),200);
  [countThis,binCentreThis] = hist(contrastThis(ii,:),200);
  count = count / sum(count);
  countThis = countThis / sum(countThis);
  aa = cumsum(count);
  interval95mass(ii,:) = [binCentre(max(find(aa < ALPHA/2))) binCentre(min(find(aa > 1 - ALPHA/2)))];
  bb = cumsum(countThis);
  interval95massThis(ii,:) = [binCentreThis(max(find(bb < ALPHA/2))) binCentreThis(min(find(bb > 1 - ALPHA/2)))];
end
subplot(3,2,4); plot(vectorTime,100*interval95mass,'-r','LineWidth',1.0); drawnow
subplot(3,2,6); plot(vectorTime,100*interval95mass,'-r','LineWidth',1.0); drawnow
subplot(3,2,3); plot(vectorTime,100*interval95massThis,'-r','LineWidth',1.0); drawnow
subplot(3,2,5); plot(vectorTime,100*interval95massThis,'-r','LineWidth',1.0); drawnow
%
%%%%%%%%%%

%%%%%%%%%%
% Function to generate convolution matrix from design matrix.
%%%%%%%%%%
function dsgMconvmtx = erConGenConvMtx(dsgM,epoch)

  dsgMconvmtx = convmtx(dsgM(:,1),epoch);
  for ii = 2:size(dsgM,2)

    dsgMconvmtx = [dsgMconvmtx convmtx(dsgM(:,ii),epoch)];
  end
  dsgMconvmtx = dsgMconvmtx(1:size(dsgM,1),:);

return;

%%%%%%%%%%
% Randomisation functions hereafter.
%%%%%%%%%%
function dsgMRand = erConRandDsg1(dsgM)
%

  dsgMRand = zeros(size(dsgM,1),size(dsgM,2));

  for ii = 1:size(dsgM,2), countsCond(ii) = sum(dsgM(:,ii)); end

  aa = randperm(size(dsgMRand,1));
  dsgMRand(aa(1:countsCond(1)),1) = 1;
  for ii = 2:size(dsgMRand,2)

    dsgMRand(aa(countsCond(ii-1)+1:countsCond(ii-1)+1+countsCond(ii)),ii) = 1;
  end

return;

function dsgMRand = erConRandDsg2(dsgM)
% This function randomises the design matrix as follows:
%   ones in the first column of the matrix stay put;
%   for all subsequent columns, ones are randomised: a one may
%   occur OFFSET_PROBE_TR TRs after the occurrence of a one in
%   the first column;
%   the number of ones in every column is preserved.

% For this regime, a probe for conditions >= 2 may only occur
% OFFSET_PROBE_TR TRs after the occurrence of a probe for
% condition 1.
OFFSET_PROBE_TR = 4;

dsgMRand = zeros(size(dsgM,1),size(dsgM,2));

for ii = 1:size(dsgM,2), countsCond(ii) = sum(dsgM(:,ii)); end

iiCandidateProbes = find(dsgM(:,1) == 1) + OFFSET_PROBE_TR;
iiCandidateProbes = iiCandidateProbes(randperm(length(iiCandidateProbes)));

dsgMRand(:,1) = dsgM(:,1);

countsCond(1) = 0; % fudge to make the following loop easy to code
for ii = 2:size(dsgMRand,2)

  dsgMRand(iiCandidateProbes((countsCond(ii-1)+1):(countsCond(ii-1)+countsCond(ii))),ii) = 1;
end

return;


