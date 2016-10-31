% modelCoranalFromGlm
%
%        $Id$
%      usage: thisView = modelCoranalFromGlm(thisView)
%         by: julien besle
%       date: 2010-10-24
%     inputs: 
%    outputs: 
% 
%    purpose: computes correlation analysis on data simulated from GLM/deconvolution estimates
%        e.g:
%
function thisView = modelCoranalFromGlm(thisView,params)


%initialize parameters
if ieNotDefined('params')
  params = modelCoranalFromGlmGUI('thisView',thisView);
  if ieNotDefined('params')
    mrWarnDlg('(modelCoranalFromGlm) Analysis cancelled');
    return;
  end
end

curGroupName = params.originalGroup;
thisView = viewSet(thisView,'curGroup',viewGet(thisView,'groupNum',curGroupName));
thisView = viewSet(thisView,'curanalysis',viewGet(thisView,'analysisNum',params.originalAnalysis));

if isempty(thisView.analyses) || ~ismember(thisView.analyses{thisView.curAnalysis}.type,{'glmAnalStats','glmAnal','erAnal','deconvAnal'})
  mrWarnDlg('You must load an event-related/GLM Analysis');
  return;
end
tic
set(viewGet(thisView,'figNum'),'Pointer','watch');drawnow;

d = viewGet(thisView,'d');
if ismember(thisView.analyses{thisView.curAnalysis}.type,{'erAnal','deconvAnal'})
  d.hrf = eye(d.hdrlen);
end
scanDims = viewGet(thisView,'dims');
curScan = viewGet(thisView,'curScan');

reverseScans = eval(params.reverseScans);
contrasts = params.contrasts;
if strcmp(params.noiseMode,'No noise')
  params.nScans = 2; %no need to model more than one sequence and its reverse
  reverseScans = 2;
  params.reverseScans = mat2str(reverseScans);
end

stimDurationTR = params.stimDurationS/d.tr;
scanDurationTR = (params.ncycles*params.cycleLengthTR);
stimCycleTR = params.cycleLengthTR/size(contrasts,1);
if strcmp(params.TRCorrection,'Yes')
  bcwdDelay = params.delayTR-1;
else
  bcwdDelay = params.delayTR;
end

nStims = size(contrasts,1);

%design a stimulus sequence with enough data on both sides to make up for dummies
totalDurationTR = scanDurationTR+2*params.dummyStimsTR;
twDesign = zeros(totalDurationTR,nStims);
%design first stim
for iCycle = 1:params.ncycles+ceil(2*params.dummyStimsTR/params.cycleLengthTR) %go a bit further to account for dummies on both sides
  twDesign((iCycle-1)*params.cycleLengthTR+(1:stimDurationTR),1)=1;
end
for iStim = 2:nStims
  twDesign(stimCycleTR*(iStim-1)+1:end,iStim) = twDesign(1:end - stimCycleTR*(iStim-1),1);
  %should we use circshift instead ? depends how the TW paradigm is designed in reality. but using dummy stims should do it
end
%cut to length + twice the dummies + any difference between stimDurationTR and stimCycleTR
twDesign = twDesign(1:totalDurationTR+stimDurationTR-stimCycleTR,1:nStims);
%make reverse sequence
reverseTwDesign = flipud(twDesign);

nHdr = size(d.ehdr,4);
nHdrComponents = size(d.ehdr,5);

co = nan(scanDims,'single');
amp = nan(scanDims,'single');
ph = nan(scanDims,'single');
averageCo = nan(scanDims,'single');
averageAmp = nan(scanDims,'single');
averagePh = nan(scanDims,'single');

%process a subset of voxels
subsetBox = eval(params.subsetBox);
subsetDims = diff(subsetBox,1,2)+1;
[maxNSlices rawNumSlices numRowsAtATime precision] = getNumSlicesAtATime(scanDurationTR*params.nScans,subsetDims,'single');
maxNSlices = min(subsetDims(3),maxNSlices);
slicesToLoad = maxNSlices*ones(1,floor(subsetDims(3)/maxNSlices));
if rem(subsetDims(3),maxNSlices)
  slicesToLoad(end+1) = rem(subsetDims(3),maxNSlices);
end
% calculate which slices we will be working on
lastSlices = subsetBox(3,1) + cumsum(slicesToLoad) -1; 
firstSlices = lastSlices - slicesToLoad +1; 

% let's not worry about junk frames for now
%junkFrames = viewGet(thisView,'totalJunkedFrames');
frameNums = [1 scanDurationTR*params.nScans];
nFrames = diff(frameNums)+1;

modelAverageTSeries = [];
for iLoad = 1:length(slicesToLoad)
  
  switch(params.noiseMode)
    case 'Residuals'
       fprintf('(loadScan) Loading scan %i from group: %s, ',curScan,curGroupName)
       tSeries = loadTSeries(thisView,curScan,[firstSlices(iLoad) lastSlices(iLoad)],frameNums,subsetBox(1,:),subsetBox(2,:), precision);
       tSeries = permute(tSeries,[4 1 2 3]);
       %convert to percent Time series
       % subtract off column means
       colmeans = repmat(mean(tSeries,1),[nFrames 1 1]);
       tSeries = tSeries - colmeans;
       % convert to percent signal change
       tSeries = 100*tSeries./colmeans;
       clear colmeans

    case 'No noise'
      fprintf('(loadScan) Processing ')
      
  end
  
  fprintf('slices=%d:%d of %d, slices=%d:%d of %d, slices=%d:%d of %d, frames=%d:%d\n',...
         firstSlices(iLoad), lastSlices(iLoad),scanDims(3),...
         subsetBox(1,1),subsetBox(1,2),scanDims(1),...
         subsetBox(2,1),subsetBox(2,2),scanDims(2),...
         frameNums(1), frameNums(2));
 
  hWaitbar = mrWaitBar(-inf,'Modeling time-series...');
  modelTSeriesTW = nan(scanDurationTR,subsetDims(1),subsetDims(2),slicesToLoad(iLoad),params.nScans,'single');
  ehdr = d.ehdr(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),firstSlices(iLoad):lastSlices(iLoad),:,:);
  %put estimates on first 2 dimensions
  ehdr = permute(ehdr,[5 4 1 2 3]);
  for z = 1:slicesToLoad(iLoad)
    fwdModelTSeries = zeros([totalDurationTR+size(d.hrf,1)-1 subsetDims(1) subsetDims(2)],'single');
    bcwdModelTSeries = zeros([totalDurationTR+size(d.hrf,1)-1 subsetDims(1) subsetDims(2)],'single');
    for y = 1:subsetDims(2)
      for x = 1:subsetDims(1)
%         if y==63 && x == 36     %DEBUG
%           keyboard
%         end
        if ~isnan(ehdr(1,1,x,y,z))

          hdr = d.hrf*ehdr(:,:,x,y,z)* contrasts';
          for iHdr = 1:nHdr
            fwdModelTSeries(:,x,y) = fwdModelTSeries(:,x,y)+convn(twDesign(:,iHdr),hdr(:,iHdr));
            bcwdModelTSeries(:,x,y) = bcwdModelTSeries(:,x,y)+convn(reverseTwDesign(:,iHdr),hdr(:,iHdr));
          end
        else
          fwdModelTSeries(:,x,y) = nan(totalDurationTR+size(d.hrf,1)-1,1,'single');
          bcwdModelTSeries(:,x,y) = nan(totalDurationTR+size(d.hrf,1)-1,1,'single');
        end
      end

    end
    %remove dummies and extra TRs from convolution
    fwdModelTSeries = fwdModelTSeries(params.dummyStimsTR+(1:scanDurationTR),:,:);
    bcwdModelTSeries = bcwdModelTSeries(params.dummyStimsTR+(1:scanDurationTR),:,:);

    fwdModelTSeries = fwdModelTSeries/nHdr;
    bcwdModelTSeries = bcwdModelTSeries/nHdr;    

    %add noise
    switch(params.noiseMode)
      case 'Residuals'
        
        tempEhdr = reshape(ehdr(:,:,:,:,z),nHdr*nHdrComponents,subsetDims(1),subsetDims(2));
        modelTSeriesGLM = nan(nFrames,subsetDims(1),subsetDims(2),'single');
        for y = 1:subsetDims(2)
          modelTSeriesGLM(:,:,y) = d.scm(frameNums(1):frameNums(2),:)*tempEhdr(:,:,y);
        end
        noise = tSeries(:,:,:,z) - modelTSeriesGLM;
        noise = permute(reshape(noise,scanDurationTR,params.nScans,subsetDims(1),subsetDims(2)),[1 3 4 2]);
        clear('modelTSeriesGLM');
        
      case 'No noise'
        noise = zeros(scanDurationTR,subsetDims(1),subsetDims(2),params.nScans,'single');
    
    end
    
    for iScan = 1:params.nScans
      if ismember(iScan,reverseScans)
        modelTSeriesTW(:,:,:,z,iScan) = bcwdModelTSeries + noise(:,:,:,iScan);
        %shift an reverse
        modelTSeriesTW(:,:,:,z,iScan) = circshift(modelTSeriesTW(:,:,:,z,iScan),bcwdDelay);
        modelTSeriesTW(:,:,:,z,iScan) = flipdim(modelTSeriesTW(:,:,:,z,iScan),1);
      else

        modelTSeriesTW(:,:,:,z,iScan) = fwdModelTSeries + noise(:,:,:,iScan);
        modelTSeriesTW(:,:,:,z,iScan) = circshift(modelTSeriesTW(:,:,:,z,iScan),params.delayTR);
      end
    end
    clear('noise');
    
    mrWaitBar(z/slicesToLoad(iLoad),hWaitbar);
  end
  clear('tSeries');
  mrCloseDlg(hWaitbar);
  
  mrDisp('Computing correlation analysis...');
  %collapse all voxels on one dimension
  modelTSeriesTW = reshape(modelTSeriesTW,scanDurationTR,subsetDims(1)*subsetDims(2)*slicesToLoad(iLoad),params.nScans);
  
  %compute the correlation analysis on the average scans
  averageModelTSeriesTW = mean(modelTSeriesTW,3);
  [coherence, amplitude, phase] = computeCoranal(averageModelTSeriesTW,params.ncycles,params.detrend,params.spatialnorm,params.trigonometricFunction);
  
  %compute the phase for each scan
  averagePhase = nan(subsetDims(1)*subsetDims(2)*slicesToLoad(iLoad),params.nScans,'single');
  for iScan = 1:params.nScans
    [dump1, dump2, averagePhase(:,iScan)] = computeCoranal(modelTSeriesTW(:,:,iScan),params.ncycles,params.detrend,params.spatialnorm,params.trigonometricFunction);
  end
  averagePhase = angle(mean(exp(1i*averagePhase),2));
  averagePhase(averagePhase<0)=averagePhase(averagePhase<0)+2*pi;
 
  %average fwd and bcwd scans separately and compute average coherence and amplitude
  bcwdModelTSeriesTW = mean(modelTSeriesTW(:,:,reverseScans),3);
  [bcwdCoherence, bcwdAmplitude] = computeCoranal(bcwdModelTSeriesTW,params.ncycles,params.detrend,params.spatialnorm,params.trigonometricFunction);
  clear('bcwdModelTSeriesTW');
  modelTSeriesTW(:,:,reverseScans) = [];
  modelTSeriesTW = mean(modelTSeriesTW,3);
  [fwdCoherence, fwdAmplitude] = computeCoranal(modelTSeriesTW,params.ncycles,params.detrend,params.spatialnorm,params.trigonometricFunction);
  clear('modelTSeriesTW')
  averageCoherence = (fwdCoherence+bcwdCoherence)/2;
  averageAmplitude = (fwdAmplitude+bcwdAmplitude)/2;

  %reshape to scan dimensions
  co(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),firstSlices(iLoad):lastSlices(iLoad)) = reshape(coherence,subsetDims(1),subsetDims(2),slicesToLoad(iLoad));
  amp(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),firstSlices(iLoad):lastSlices(iLoad)) = reshape(amplitude,subsetDims(1),subsetDims(2),slicesToLoad(iLoad));
  ph(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),firstSlices(iLoad):lastSlices(iLoad)) = reshape(phase,subsetDims(1),subsetDims(2),slicesToLoad(iLoad));
  averageCo(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),firstSlices(iLoad):lastSlices(iLoad)) = reshape(averageCoherence,subsetDims(1),subsetDims(2),slicesToLoad(iLoad));
  averageAmp(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),firstSlices(iLoad):lastSlices(iLoad)) = reshape(averageAmplitude,subsetDims(1),subsetDims(2),slicesToLoad(iLoad));
  averagePh(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),firstSlices(iLoad):lastSlices(iLoad)) = reshape(averagePhase,subsetDims(1),subsetDims(2),slicesToLoad(iLoad));
  modelAverageTSeries = cat(4,modelAverageTSeries,reshape(averageModelTSeriesTW,scanDurationTR,subsetDims(1),subsetDims(2),slicesToLoad(iLoad)));
  
  fprintf('\n');

end

%Create Group
if isempty(viewGet(thisView,'groupnum',params.newGroupName))
  thisView = viewSet(thisView,'newgroup',params.newGroupName);
end


% info for new saved scan (the modeled average Tseries
%(using header of current scan as the template for the nifti header)
scanParams.originalFileName{1} = viewGet(thisView,'tSeriesFile',curScan);
scanParams.originalGroupName{1} = curGroupName;
scanParams.fileName = [];
scanParams.junkFrames = 0;
scanParams.nFrames = scanDurationTR;
scanParams.description = ['average modeled from ' params.originalAnalysis ...
  ' - ' num2str(params.nScans) ' scans (' mat2str(reverseScans) ' reversed) - delay=' num2str(params.delayTR)];
if strcmp(params.TRCorrection,'Yes')
  scanParams.description = [scanParams.description ' - TR corrected'];
end
if strcmp(params.noiseMode,'No noise')
  scanParams.description = [scanParams.description ' - No Noise'];
end
scanParams.vol2mag = viewGet(thisView,'scanVol2mag',curScan);
scanParams.vol2tal = viewGet(thisView,'scanVol2tal',curScan);
hdr = cbiReadNiftiHeader(viewGet(thisView,'tseriesPath',curScan));


%the following has to be done in this specific order:

%first set the group where we're gonna put everything
newGroupNum = viewGet(thisView,'groupNum',params.newGroupName);
thisView = viewSet(thisView,'currentGroup',newGroupNum);

%load the analysis and overlays or create new ones (with old number of scans)
dateString = datestr(now);
corAnal = viewGet(thisView,'analysis',viewGet(thisView,'analysisNum',params.saveName));
%if there is an analysis with this name but it hasn't been created with modelCoranalFromGlm, change the name
if ~isempty(corAnal) && ~strcmp(corAnal.function,'modelCoranalFromGlm')
  mrWarnDlg('An analysis with this name already exists but had not been created with modelCoranalFromGlm');
  corAnal = [];
  params.saveName = [params.saveName '2'];
end
if isempty(corAnal) 
  nScans = viewGet(thisView,'nScans');
  % params for scans with no analysis/overlays
  emptyParams.scanList = 1:nScans;
  emptyParams.scanParams = cell(1,nScans);
  for iScan = 1:nScans 
    emptyParams.tseriesFile{iScan} = viewGet(thisView,'tseriesFile',iScan,newGroupNum);
  end
  %Create Analysis
  corAnal.name = params.saveName; 
  corAnal.type = 'corAnal';
  corAnal.groupName = params.newGroupName;
  corAnal.function = 'modelCoranalFromGlm';
  corAnal.reconcileFunction = 'defaultReconcileParams';
  corAnal.mergeFunction = 'defaultMergeParams';
  corAnal.guiFunction = 'modelCoranalFromGlmGUI';
  corAnal.params = emptyParams;
  %corAnal.overlayInterpFunction = 'corAnalInterp'; this does not exist
  corAnal.date = dateString;
  thisView = viewSet(thisView,'newanalysis',corAnal);
  
  %create overlay template
  coOverlay.groupName = params.newGroupName;
  coOverlay.function = 'modelCoranalFromGlm';
  coOverlay.reconcileFunction = 'defaultReconcileParams';
  coOverlay.mergeFunction = 'defaultMergeParams';
  coOverlay.type = 'corAnal';
  coOverlay.interrogator = 'corAnalPlot';
  coOverlay.date = dateString;
  coOverlay.params = emptyParams;
  coOverlay.data = cell(1,nScans);
  
  %coherence
  coOverlay.colormap = jet(256);
  coOverlay.range = [0 1];
  coOverlay.clip = [0 1];
  coOverlay.name = 'co';
  %average coherence
  averageCoOverlay = coOverlay;
  averageCoOverlay.name = 'co (average)';
  %amplitude
  ampOverlay=coOverlay;
  ampOverlay.name = 'amp';
  ampOverlay.colormap = hot(256);
  ampOverlay.range = [0 10];
  ampOverlay.clip = [0 inf];
  %average amplitude
  averageAmpOverlay=ampOverlay;
  averageAmpOverlay.name = 'amp (average)';
  %phase
  phOverlay = coOverlay;
  phOverlay.name = 'ph';
  phOverlay.range = [0 2*pi];
  phOverlay.clip = [0 2*pi];
  phOverlay.colormap = hsv(256);
  %average phase
  averagePhOverlay=phOverlay;
  averagePhOverlay.name = 'ph (average)';

else
  thisView = viewSet(thisView,'curanalysis',viewGet(thisView,'analysisnum',params.saveName));
  coOverlay = viewGet(thisView,'overlay',viewGet(thisView,'overlaynum','co'));
  phOverlay = viewGet(thisView,'overlay',viewGet(thisView,'overlaynum','ph'));
  ampOverlay = viewGet(thisView,'overlay',viewGet(thisView,'overlaynum','amp'));
  averageAmpOverlay = viewGet(thisView,'overlay',viewGet(thisView,'overlaynum','amp (average)'));
  averageCoOverlay = viewGet(thisView,'overlay',viewGet(thisView,'overlaynum','co (average)'));
  averagePhOverlay = viewGet(thisView,'overlay',viewGet(thisView,'overlaynum','ph (average)'));
end


%then add the new scan
%fill the new scan with data
newTSeries = nan(scanDims(1),scanDims(2),scanDims(3),scanDurationTR,'single');
newTSeries(subsetBox(1,1):subsetBox(1,2),subsetBox(2,1):subsetBox(2,2),subsetBox(3,1):subsetBox(3,2),:) = permute(modelAverageTSeries,[2 3 4 1]);
clear modelAverageTSeries
thisView = saveNewTSeries(thisView,newTSeries,scanParams,hdr);
clear('newTSeries')
lastScan = viewGet(thisView,'nScans');
tseriesFileName = viewGet(thisView,'tseriesfile',lastScan);

%and add the overlays (with the new number of scans)
%this can only be done once the number of scans matches the new number of overlays
coOverlay.data{lastScan} = co;
coOverlay.params.scanParams{lastScan} = params;
coOverlay.params.scanList(lastScan) = lastScan;
coOverlay.params.tseriesFile{lastScan} = tseriesFileName;
thisView = viewSet(thisView,'newoverlay',coOverlay);  
ampOverlay.data{lastScan} = amp;
ampOverlay.params.scanParams{lastScan} = params;
ampOverlay.params.scanList(lastScan) = lastScan;
ampOverlay.params.tseriesFile{lastScan} = tseriesFileName;
thisView = viewSet(thisView,'newoverlay',ampOverlay);  
phOverlay.data{lastScan} = ph;
phOverlay.params.scanParams{lastScan} = params;
phOverlay.params.scanList(lastScan) = lastScan;
phOverlay.params.tseriesFile{lastScan} = tseriesFileName;
thisView = viewSet(thisView,'newoverlay',phOverlay);  
averageCoOverlay.data{lastScan} = averageCo;
averageCoOverlay.params.scanParams{lastScan} = params;
averageCoOverlay.params.scanList(lastScan) = lastScan;
averageCoOverlay.params.tseriesFile{lastScan} = tseriesFileName;
thisView = viewSet(thisView,'newoverlay',averageCoOverlay);  
averageAmpOverlay.data{lastScan} = averageAmp;
averageAmpOverlay.params.scanParams{lastScan} = params;
averageAmpOverlay.params.scanList(lastScan) = lastScan;
averageAmpOverlay.params.tseriesFile{lastScan} = tseriesFileName;
thisView = viewSet(thisView,'newoverlay',averageAmpOverlay);  
averagePhOverlay.data{lastScan} = averagePh;
averagePhOverlay.params.scanParams{lastScan} = params;
averagePhOverlay.params.scanList(lastScan) = lastScan;
averagePhOverlay.params.tseriesFile{lastScan} = tseriesFileName;
thisView = viewSet(thisView,'newoverlay',averagePhOverlay);  

%change analysis parameters to match the number of overlays
%this can only be done once the number of scans and overlays match the new number of parameters
corAnal = viewGet(thisView,'analysis',viewGet(thisView,'analysisNum',params.saveName));
corAnal.params.tseriesFile{lastScan} = tseriesFileName;
corAnal.params.scanParams{lastScan} = params;
corAnal.params.scanList(lastScan) = lastScan;
thisView = viewSet(thisView,'newanalysis',corAnal);

saveAnalysis(thisView,corAnal.name);

set(viewGet(thisView,'figNum'),'Pointer','arrow');drawnow

refreshMLRDisplay(viewGet(thisView,'viewNum'));

toc

