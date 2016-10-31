% glmAnalysis.m
%
%      usage: thisView = glmAnalysis(thisView,params)
%         by: farshad moradi,  modified by julien besle 12/01/2010 to 25/10/2010 to perform statistical inference on parameter/contrast estimates
%       date: 06/14/07
%    purpose: GLM analysis using design matrix convolved with any HRF model (including deconvolution)
%              $Id$
%
% 
%     Fits a GLM to the data using ordinary or generalized Least Squares
%     GLM is composed of EVs which can be stimulus times or combinations thereof
%     outputs r2 overlay
%     computes any number of contrasts
%     T-tests are performed on any linear combinations of Explanatory Variables (EVs) against zero (H0: contrast'*beta==0)
%     F-tests test overall effect of any set of contrasts, not necessarily identical to above contrasts (H0: no contrast differs from zero)

function [thisView,params] = glmAnalysis(thisView,params,varargin)

% check arguments
if ~any(nargin == [1 2 3 4 5])
  help glmAnalysis
  return
end

mrGlobals;

% other arguments
eval(evalargs(varargin));
if ieNotDefined('justGetParams'),justGetParams = 0;end
if ieNotDefined('defaultParams'),defaultParams = 0;end
if ieNotDefined('scanList'),scanList = [];end
if ieNotDefined('params'),params = [];end

% First get parameters
if isempty(params) || justGetParams
  params = glmAnalysisGUI('thisView',thisView,'params',params,'defaultParams',defaultParams,'scanList',scanList);
end

% Abort if params empty
if ieNotDefined('params')
  disp('(glmAnalysis) GLM analysis cancelled');
  return
% just return parameters
elseif justGetParams
  return
end


% set the group
thisView = viewSet(thisView,'groupName',params.groupName);
% Reconcile params with current status of group and ensure that it has
% the required fields. 
params = defaultReconcileParams([],params);  %this seems quite useless to me
%is it just here to link parameters to file names ?

%just to have shorter variable names
scanParams = params.scanParams;
numberFtests = length(params.restrictions);
numberContrasts = size(params.contrasts,1);
numberTests = numberFtests+numberContrasts;
computePermutations = numberTests && (params.permutationTests || (params.parametricTests && params.permutationFweAdjustment));

if params.covCorrection   %number of voxels to get around the ROI/subset box in case the covariance matrix is estimated
  voxelsMargin = repmat(floor(params.covEstimationAreaSize/2),1,3);
  switch(params.covEstimationPlane)
    case {'Sagittal'}
      voxelsMargin(1)=0;
    case {'Axial'}
      voxelsMargin(3)=0;
    case {'Coronal'}
      voxelsMargin(2)=0;
  end
else
   voxelsMargin = [0 0 0];
end
if params.spatialSmoothing  %we'll also need a margin if we're spatially smoothing
  switch(params.smoothingPlane)
      case {'Sagittal'}
        voxelsMargin = max(voxelsMargin, [0 params.spatialSmoothing params.spatialSmoothing]);
      case {'Axial'}
        voxelsMargin = max(voxelsMargin,[params.spatialSmoothing params.spatialSmoothing 0]);
      case {'Coronal'}
        voxelsMargin = max(voxelsMargin,[params.spatialSmoothing 0 params.spatialSmoothing]);
      case '3D'
        voxelsMargin = max(voxelsMargin,repmat(params.spatialSmoothing,1,3));
  end
end
%--------------------------------------------------------- Main loop over scans ---------------------------------------------------
set(viewGet(thisView,'figNum'),'Pointer','watch');drawnow;
%initialize the data we're keeping for output overlays
precision = mrGetPref('defaultPrecision');
r2 = cell(1,params.scanNum(end));
if numberContrasts
  contrast = r2;
end
if numberTests
  if params.parametricTests
    parametricP = r2;  
    if params.outputStatistic
      statistic = r2;
    end
    fdrParametricP=r2;
    fweParametricP=r2;
    if params.bootstrapFweAdjustment
      bootstrapFweParametricP = r2;
    end
    if params.permutationFweAdjustment
      permuteFweParametricP = r2;
    end
    if params.TFCE 
      if params.bootstrapFweAdjustment
        bootstrapFweTfceP = r2;
      end
      if params.permutationFweAdjustment
        permuteFweTfceP = r2;
      end
    end
  end
  if params.permutationTests
    permuteP = r2;
    fdrPermuteP=r2;
    fwePermuteP=r2;
    if params.TFCE
      tfcePermuteP = r2;
      fdrTfcePermuteP=r2;
      fweTfceRandT=r2;
    end
  end
  if params.bootstrapTests
    bootstrapP = r2;
    fdrBootstrapP=r2;
    fweBootstrapP=r2;
    if params.bootstrapFweAdjustment && ~params.noDoubleBootstrap
      bootstrapFweBootstrapP = r2;
    end
    if params.TFCE  
      tfceBootstrapP = r2;
      fdrTfceBootstrapP=r2;
      fweTfceBootstrapP=r2;
      if params.bootstrapFweAdjustment && ~params.noDoubleBootstrap
        bootstrapFweTfceBootstrapP = r2;
      end
    end
  end
end


for iScan = params.scanNum
  numVolumes = viewGet(thisView,'nFrames',iScan);
  scanDims{iScan} = viewGet(thisView,'dims',iScan);

  %compute the dimensions of the subset of voxels to load
  switch(params.analysisVolume)
   case {'Whole volume','Subset box'}
       subsetBox{iScan} = eval(scanParams{iScan}.subsetBox);
   case {'Loaded ROI(s)','Visible ROI(s)'}
      %get the smallest box containing all the voxels from all the (visible) ROIs 
      if strcmp(params.analysisVolume,'Visible ROI(s)')
        roiList = viewGet(thisView,'visibleRois');
      else
        roiList = 1:viewGet(thisView,'numberOfRois');
      end
      [subsetBox{iScan}, whichRoi, marginVoxels] = getRoisBox(thisView,iScan,voxelsMargin,roiList);
      usedVoxelsInBox = marginVoxels | any(whichRoi,4);
      %clear('whichRoi','marginVoxels');
      if params.covCorrection && ~strcmp(params.covEstimationBrainMask,'None')
        [dump,brainMaskRoiNum] = ismember(params.covEstimationBrainMask,viewGet(thisView,'roiNames'));
        brainMaskScanCoords = getROICoordinates(thisView,brainMaskRoiNum,iScan);
        %keep only those voxels that are in the subset box
        brainMaskScanCoords(:,any(brainMaskScanCoords-repmat(subsetBox{iScan}(:,2),1,size(brainMaskScanCoords,2))>0))=[];
        brainMaskScanCoords=brainMaskScanCoords-repmat(subsetBox{iScan}(:,1),1,size(brainMaskScanCoords,2))+1;
        brainMaskScanCoords(:,any(brainMaskScanCoords<1))=[];
        if ~isempty(brainMaskScanCoords)
          %create a mask the same size as the subsetbox
          covEstimationBrainMask = false(size(usedVoxelsInBox));
          covEstimationBrainMask(sub2ind(size(usedVoxelsInBox),brainMaskScanCoords(1,:)',brainMaskScanCoords(2,:)',brainMaskScanCoords(3,:)'))=true;
        else
          covEstimationBrainMask = [];
        end
      end
  end
  subsetDims = diff(subsetBox{iScan},1,2)'+1;
  
  %Compute the number of slices to load at once in order to minimize memory usage
  if params.TFCE && computePermutations
    %in this particular case, we force loading the whole dataset at once
    slicesToLoad = subsetDims(3);
    loadCallsPerBatch = {1};
    rawNumSlices = subsetDims(3);
    %compare the size of the array to load to the memory preference and ask for confirmation if larger than maxBlockSize
    switch (precision)
      case {'double'}
        bytesPerNum = 8;
      case{'single'}
        bytesPerNum = 4;
    end
    sizeArray = bytesPerNum*numVolumes*prod(subsetDims);
    if sizeArray> mrGetPref('maxBlocksize')
      if ~askuser(sprintf('(glmAnalysis) This will load an array of %.2f Gb in memory. Are you sure you want to proceed ?', sizeArray/1024^3));
        return;
      end
    end
  else
    % otherwise, see how many slices we can load at once 
    [maxNSlices rawNumSlices numRowsAtATime precision] = getNumSlicesAtATime(numVolumes,subsetDims,precision);
    maxNSlices = min(subsetDims(3),maxNSlices);
    slicesToLoad = maxNSlices*ones(1,floor(subsetDims(3)/maxNSlices));
    if rem(subsetDims(3),maxNSlices)
      slicesToLoad(end+1) = rem(subsetDims(3),maxNSlices);  %list of the number of slices to load at once
    end
    
    switch(params.analysisVolume)
      case {'Loaded ROI(s)','Visible ROI(s)'}
        % for ROIs, since the box contains voxels we won't use, 
        % choose how many loads we can do while still computing all voxels at once (=batch)
        %(we assume that we can at least load one entire slice by calls to loadScan)
        loadCallsPerBatch = {[]};   %this will contains nBatches lists of loads (indices to the slicesToLoad vector) whose total number of useful voxels can be held in memory at once
        nVoxels = 0;
        currentFirstSlice = 0;
        nBatches = 1;
        for iLoad = 1:length(slicesToLoad)
          thisNvoxels = nnz(usedVoxelsInBox(:,:,currentFirstSlice+(1:slicesToLoad(iLoad))) > 0);
          nVoxels =  nVoxels+thisNvoxels; %total number of useful voxels loaded in this batch
          % check if the total number of voxels can be held at once in memory
          [roiSlicesAtATime rawRoiSlices] = getNumSlicesAtATime(numVolumes,[nVoxels 1 1],precision);
          if rawRoiSlices<1 %if the returned value is less than 1, this means it can't, and we have to create another batch
            nBatches = nBatches+1;
            currentFirstSlice = sum(slicesToLoad(1:iLoad));
            loadCallsPerBatch{nBatches} = [];
            nVoxels=thisNvoxels;
          end
          loadCallsPerBatch{nBatches} = [loadCallsPerBatch{nBatches} iLoad];  %add the current load to the list of loads in the batch
        end
        if nBatches>1 && ( (params.spatialSmoothing && ~strcmp(params.smoothingPlane,'Axial')) || (params.covCorrection && ~strcmp(params.covEstimationPlane,'Axial')))
          mrWarnDlg('(glmAnalysis) Smoothing or covariance estimation planes other than ''Axial'' cannot be used if the analysis is run separately on subsets of slices because of memory limits.');
          mrWarnDlg('Increase memory block size or reduce the size of the ROI(s). Aborting...');
          return;
        end

      case {'Whole volume','Subset box'}
        loadCallsPerBatch = num2cell(1:length(slicesToLoad)); %in this case, each batch comprises of only one load

    end
  end

  if rawNumSlices<1
      mrWarnDlg(['Too many data points to perform the analysis on whole slices. Implement row analysis support, increase memory block size or reduce subset box/ROI(s) size (by a factor ' num2str(subsetDims(3)/rawNumSlices) '). Aborting...']);
      return;
  end

  ehdr = []; s2 = []; autoCorrelationParameters = [];
  ehdrBootstrapCIs = [];  contrastBootstrapCIs = [];


  % calculate which slices we will be working on
  lastSlices = cumsum(slicesToLoad); 
  firstSlices = lastSlices - slicesToLoad + 1; 
  %loop
  for iBatch = 1:length(loadCallsPerBatch)
    switch(params.analysisVolume)
      case {'Whole volume','Subset box'}
        %Load data
        d = loadScan(thisView,iScan,[],subsetBox{iScan}(3,1) + [firstSlices(iBatch) lastSlices(iBatch)] -1, precision,subsetBox{iScan}(1,:),subsetBox{iScan}(2,:));
      case {'Loaded ROI(s)','Visible ROI(s)'}
        for iLoad= loadCallsPerBatch{iBatch}
          dummy = loadScan(thisView,iScan,[],subsetBox{iScan}(3,1) + [firstSlices(iLoad) lastSlices(iLoad)] -1, precision,subsetBox{iScan}(1,:),subsetBox{iScan}(2,:));
          dummy.data = reshape(dummy.data,[prod(dummy.dim(1:3)) dummy.dim(4)]);
          dummy.data = dummy.data(usedVoxelsInBox(:,:,firstSlices(iLoad):lastSlices(iLoad))>0,:);
          if iLoad==loadCallsPerBatch{iBatch}(1)
            d=dummy;
          else
            d.data = cat(1,d.data, dummy.data);
          end
        end
        clear('dummy');
        d.data = permute(d.data,[1 3 4 2]);
    end
    clear('dummy');
    
    if iBatch==1
      %----------------------------------------Design matrix: same for all voxels/runs
      % get the stim volumes, if empty then abort
      d = getStimvol(d,scanParams{iScan});
      if isempty(d.stimvol),mrWarnDlg('No stim volumes found');return,end
      % do any call for preprocessing
      if ~isempty(scanParams{iScan}.preprocess)
        d = eventRelatedPreProcess(d,scanParams{iScan}.preprocess);
      end

      actualStimvol = d.stimvol;
      %precompute permutation vectors
      if computePermutations
        nResamples = params.nResamples;
        %find which event types will be randomized
        %first which EVs are involved in all contrasts/restrictions
        contrastEVs = any([params.contrasts;cell2mat(params.restrictions)],1);
        %then which event types constitute these EVs
        randEvents = any(scanParams{iScan}.stimToEVmatrix(:,contrastEVs),2);
        %now compute permutation indices for those events, while keeping other events indices unchanged
        nEventTypes = length(actualStimvol);
        numberEvents = NaN(nEventTypes,1);
        totalNumberEvents = 0;
        stimvolToPermuteIndices = [];
        for iEventType = 1:nEventTypes
           numberEvents(iEventType) = length(actualStimvol{iEventType});
           if randEvents(iEventType)
             stimvolToPermuteIndices = [stimvolToPermuteIndices totalNumberEvents+(1:numberEvents(iEventType))];
           end
           totalNumberEvents = totalNumberEvents+numberEvents(iEventType);
        end
        stimvolToPermute = cell2mat(actualStimvol(randEvents));
        nStimsToPermute = length(stimvolToPermute);
        permutations = repmat(cell2mat(actualStimvol),nResamples,1);
        for iPerm = 1:nResamples
           permutations(iPerm,stimvolToPermuteIndices) = stimvolToPermute(randperm(nStimsToPermute));
        end
      else
        nResamples = 0;
      end

      %create model HRF
      [params.hrfParams,d.hrf] = feval(params.hrfModel, params.hrfParams, d.tr/d.designSupersampling,scanParams{iScan}.acquisitionDelay,1);

      d.volumes = 1:d.dim(4);
      %make a copy of d
      actualD = d;

    else
      thisD=d;
      d = actualD;
      d.data = thisD.data;
      clear('thisD');
    end
    d.dim = size(d.data);
    
    
    switch(params.analysisVolume)
      case {'Loaded ROI(s)','Visible ROI(s)'}
            d.roiPositionInBox = any(whichRoi(:,:,firstSlices(loadCallsPerBatch{iBatch}(1)):lastSlices(loadCallsPerBatch{iBatch}(end)),:),4);
        d.marginVoxels = marginVoxels(:,:,firstSlices(loadCallsPerBatch{iBatch}(1)):lastSlices(loadCallsPerBatch{iBatch}(end)),:);
        if params.covCorrection && ~strcmp(params.covEstimationBrainMask,'None')&& ~isempty(covEstimationBrainMask)
          d.covEstimationBrainMask = covEstimationBrainMask(:,:,firstSlices(loadCallsPerBatch{iBatch}(1)):lastSlices(loadCallsPerBatch{iBatch}(end)),:);
        else
          d.covEstimationBrainMask = [];
        end
    end
      
    %-------------------------------Permutations------------------------------------------
    for iPerm = 1:nResamples+1
      if iPerm==2
        hWaitBar = mrWaitBar(0,['(glmAnalysis) Computing Permutations for scan ' num2str(iScan)]);
      end
      if iPerm == 1 %if it is the actual data
        d.stimvol = actualStimvol;
        actualData = 1;
        verbose = 1;
      else
        actualData = 0;
        mrWaitBar(iPerm/nResamples,hWaitBar);
        permutedStimvol = permutations(iPerm-1,:);
        for stimnum = 1:length(actualStimvol)%randomize stim events
           d.stimvol{stimnum} = permutedStimvol(1:numberEvents(stimnum));
           permutedStimvol = permutedStimvol(numberEvents(stimnum)+1:end);
        end
        verbose = 0;
      end

      %compute the design matrix for this permutation
      d = makeDesignMatrix(d,params,verbose, iScan);
      d.emptyEVcomponents = testDesignMatrix(d.scm,d.nhdr,d.nHrfComponents,params.EVnames);
      if ~isempty(d.emptyEVcomponents)
        % for deconvolution with design supersampling, some components might be impossible to estimate
        mrWarnDlg(sprintf('(glmAnalysis) Not enough data in scan %i to estimate all EV components, ignoring empty components (This has not been tested...)',iScan));
        %they will be removed from the design matrix in getGlmStatistics and their parameter estimates replaced by NaNs
      end

      % compute estimates and statistics
      [d, out] = getGlmStatistics(d, params, verbose, precision, actualData);%, computeTtests,computeBootstrap);
      
    
      if iPerm==1
        %keep values for permutation tests
        if computePermutations
          actualStatistic = out.statistic;
          if params.permutationTests
            permuteStatisticCount = zeros(size(actualStatistic),precision); %will count how many permutation values are above the actual statistic value 
            permuteStatisticCount(isnan(actualStatistic)) = NaN; 
          end
          if params.parametricTests && params.permutationFweAdjustment
            permuteFweStatisticCount = zeros(size(actualStatistic),precision); %same for the max permutation value across space, for FWE adjustment
            permuteFweStatisticCount(isnan(actualStatistic)) = NaN; 
            permuteFweData = initResampleFWE(reshape(actualStatistic,prod(d.dim(1:3)),numberTests),params,out.parametricP);
          end
        end
        
        %reshape data into the subsetbox
        if ismember(params.analysisVolume,{'Loaded ROI(s)' 'Visible ROI(s)'}) 
          d.ehdr = reshapeToRoiBox(d.ehdr,d.roiPositionInBox);
          out.r2 = reshapeToRoiBox(out.r2,d.roiPositionInBox);
          d.s2 = reshapeToRoiBox(d.s2,d.roiPositionInBox);
          if params.covCorrection
            d.autoCorrelationParameters = reshapeToRoiBox(d.autoCorrelationParameters,d.roiPositionInBox);
          end
          if params.bootstrapTests && params.bootstrapIntervals
            d.ehdrBootstrapCIs = reshapeToRoiBox(d.ehdrBootstrapCIs,d.roiPositionInBox);
          end
          if numberTests
            if numberContrasts && (nnz(params.componentsToTest)==1 || strcmp(params.componentsCombination,'Add'))
              out.contrast = reshapeToRoiBox(out.contrast,d.roiPositionInBox);
              if params.bootstrapTests && params.bootstrapIntervals
                d.contrastBootstrapCIs = reshapeToRoiBox(d.contrastBootstrapCIs,d.roiPositionInBox);
              end
            end
            if params.parametricTests
              out.parametricP = reshapeToRoiBox(out.parametricP,d.roiPositionInBox);
              if params.outputStatistic
                out.statistic = reshapeToRoiBox(out.statistic,d.roiPositionInBox);
              end
              if params.bootstrapFweAdjustment
                out.bootstrapFweParametricP = reshapeToRoiBox(out.bootstrapFweParametricP,d.roiPositionInBox);
                if params.TFCE 
                  out.bootstrapFweTfceP = reshapeToRoiBox(out.bootstrapFweTfceP,d.roiPositionInBox);
                end
              end
            end
            if params.bootstrapTests
              out.bootstrapP = reshapeToRoiBox(out.bootstrapP,d.roiPositionInBox);
              if params.bootstrapFweAdjustment && ~params.noDoubleBootstrap
                out.bootstrapFweBootstrapP = reshapeToRoiBox(out.bootstrapFweBootstrapP,d.roiPositionInBox);
              end
              if params.TFCE 
                out.tfceBootstrapP = reshapeToRoiBox(out.tfceBootstrapP,d.roiPositionInBox);
                if params.bootstrapFweAdjustment && ~params.noDoubleBootstrap
                   out.bootstrapFweTfceBootstrapP = reshapeToRoiBox(out.bootstrapFweTfceBootstrapP,d.roiPositionInBox);
                end
              end
            end
          end
        end
        
        %concatenate data
        ehdr = cat(3,ehdr,d.ehdr);
        r2{iScan} = cat(3,r2{iScan},out.r2);
        s2 = cat(3,s2,d.s2);
        if params.covCorrection
          autoCorrelationParameters = cat(3,autoCorrelationParameters,d.autoCorrelationParameters);
        end
        if params.bootstrapTests && params.bootstrapIntervals
          ehdrBootstrapCIs = cat(3,ehdrBootstrapCIs,d.ehdrBootstrapCIs);
        end
        if numberTests
          if numberContrasts
            if nnz(params.componentsToTest)==1 || strcmp(params.componentsCombination,'Add')
              contrast{iScan} = cat(3,contrast{iScan},out.contrast);
            end
            if params.bootstrapTests && params.bootstrapIntervals  && (nnz(params.componentsToTest)==1 || strcmp(params.componentsCombination,'Add'))
              contrastBootstrapCIs = cat(3,contrastBootstrapCIs,d.contrastBootstrapCIs);
            end
          end
          if (numberFtests ||params.computeTtests) && params.parametricTests
            parametricP{iScan} = cat(3,parametricP{iScan},out.parametricP);
            if params.outputStatistic
              statistic{iScan} = cat(3,statistic{iScan},out.statistic);
            end
            if params.bootstrapFweAdjustment
              bootstrapFweParametricP{iScan} = cat(3,bootstrapFweParametricP{iScan},out.bootstrapFweParametricP);
              if params.TFCE 
                bootstrapFweTfceP{iScan} = cat(3,bootstrapFweTfceP{iScan},out.bootstrapFweTfceP);
              end
            end
          end
          if params.bootstrapTests
            bootstrapP{iScan} = cat(3,bootstrapP{iScan},out.bootstrapP);
            if params.bootstrapFweAdjustment && ~params.noDoubleBootstrap
              bootstrapFweBootstrapP{iScan} = cat(3,bootstrapFweBootstrapP{iScan},out.bootstrapFweBootstrapP);
            end
            if params.TFCE 
              tfceBootstrapP{iScan} = cat(3,tfceBootstrapP{iScan},out.tfceBootstrapP);
              if params.bootstrapFweAdjustment && ~params.noDoubleBootstrap
                bootstrapFweTfceBootstrapP{iScan} = cat(3,bootstrapFweTfceBootstrapP{iScan},out.bootstrapFweTfceBootstrapP);
              end
            end
          end
        end
        %We compute the TFCE here only in the case we forced loading the whole volume at once
        if params.TFCE && computePermutations
          thisTfce = applyFslTFCE(out.statistic,'',0);  
          %put NaNs back
          thisTfce(isnan(out.statistic)) = NaN;
          actualTfce = thisTfce;           %keep the TFCE transform
          if params.permutationTests
            permuteTfceCount = zeros(size(out.statistic),precision);     %these will count how many permutation values are above the actual TFCE values voxelwise
            permuteTfceCount(isnan(out.statistic)) = NaN;%put NaNs back
          end
          if params.parametricTests && params.permutationFweAdjustment
            permuteFweTfceCount = zeros(size(out.statistic),precision); %same for permutation adjusted FWE
            permuteFweTfceCount(isnan(out.statistic)) = NaN;%put NaNs back
            permuteFweTfceData = initResampleFWE(reshape(actualTfce,numel(thisTfce(:,:,:,1)),numberContrasts+numberFtests),params);
          end
          clear('out');
        end
      
      else % iPerm>1
        if params.permutationTests
          permuteStatisticCount = permuteStatisticCount+double(out.statistic>actualStatistic);
        end
        if params.parametricTests && params.permutationFweAdjustment
          permuteFweStatisticCount = permuteFweStatisticCount + reshape(resampleFWE(reshape(out.statistic,prod(d.dim(1:3)),numberTests),...
                                                            reshape(actualStatistic,prod(d.dim(1:3)),numberTests),...
                                                            permuteFweData), [d.dim(1:3) numberTests]);
        end
        if params.TFCE
          if isfield(d,'roiPositionInBox') 
            out.statistic = reshapeToRoiBox(out.statistic,d.roiPositionInBox); 
          end
          %We compute the TFCE here only in the case we forced loading the whole volume at once
          %NaNs in the data will be transformed to 0 by FSL, so for ROIs, 
          %TFCE is applied to the smallest box including all the ROIs, but replacing non-computed data by 0
          thisTfce = applyFslTFCE(out.statistic,'',0);  
          %put NaNs back
          thisTfce(isnan(out.statistic)) = NaN;
          clear('out');
          if params.permutationTests
            permuteTfceCount = permuteTfceCount+ double(thisTfce>actualTfce);
          end
          if params.parametricTests && params.permutationFweAdjustment
            %permuteFweTfceCount = permuteFweTfceCount+ double(repmat(max(max(max(thisTfce,[],3),[],2),[],1),[subsetDims 1])>actualTfce);
            permuteFweTfceCount = permuteFweTfceCount + reshape(resampleFWE(reshape(thisTfce,numel(thisTfce(:,:,:,1)),numberContrasts+numberFtests),...
                                                                reshape(actualTfce,numel(thisTfce(:,:,:,1)),numberContrasts+numberFtests),...
                                                                permuteFweTfceData),size(permuteFweTfceCount));
          end
          clear('thisTfce');
        end
      end
    end
    if computePermutations
      clear('actualStatistic','actualTfce')
      mrCloseDlg(hWaitBar);
    end
    d = rmfield(d,'data');

    if computePermutations %compute P-values, reshape and concatenate if needed
      if params.permutationTests
        permuteStatisticCount = computeRandP(permuteStatisticCount,nResamples);
        if ismember(params.analysisVolume,{'Loaded ROI(s)' 'Visible ROI(s)'}) %reshape data into the subsetbox
          permuteStatisticCount = reshapeToRoiBox(permuteStatisticCount,d.roiPositionInBox);
        end
        permuteP{iScan} = cat(3,permuteP{iScan},permuteStatisticCount);  
        clear('permuteStatisticCount');
        if params.TFCE 
          tfcePermuteP{iScan} = computeRandP(permuteTfceCount,nResamples);
          clear('permuteTfceCount');
          tfcePermuteP{iScan}(isnan(permuteP{iScan})) = NaN; %put NaNs back in place
        end
      end
      if params.parametricTests && params.permutationFweAdjustment
        permuteFweStatisticCount = computeRandP(permuteFweStatisticCount,nResamples);
        permuteFweStatisticCount(isnan(permuteFweStatisticCount)) = NaN; %put NaNs back in place
        if ismember(params.analysisVolume,{'Loaded ROI(s)' 'Visible ROI(s)'}) %reshape data into the subsetbox
          permuteFweStatisticCount = reshapeToRoiBox(permuteFweStatisticCount,d.roiPositionInBox);
        end
        permuteFweParametricP{iScan} = cat(3,permuteFweParametricP{iScan},permuteFweStatisticCount);
        clear('permuteFweStatisticCount');
        if params.TFCE 
          permuteFweTfceP{iScan} = computeRandP(permuteFweTfceCount,nResamples);
          clear('permuteFweTfceCount');
          permuteFweTfceP{iScan} = reshape(enforceMonotonicityResampleFWE(...
                                      reshape(permuteFweTfceP{iScan},numel(permuteFweTfceP{iScan}(:,:,:,1)),numberContrasts+numberFtests),...
                                      permuteFweTfceData), size(permuteFweTfceP{iScan}));
          permuteFweTfceP{iScan}(isnan(permuteP{iScan})) = NaN; %put NaNs back in place
        end
      end
      d = rmfield(d,'roiPositionInBox');
    end
  end
  clear ('usedVoxelsInBox');
  
  if (numberFtests ||params.computeTtests) && numberTests
    %make parametric probability maps
    if params.parametricTests
      [parametricP{iScan},fdrParametricP{iScan},fweParametricP{iScan}] = transformStatistic(parametricP{iScan},precision,params); 
      if params.permutationFweAdjustment
        permuteFweParametricP{iScan} = transformStatistic(permuteFweParametricP{iScan},precision,params); 
        if params.TFCE 
          permuteFweTfceP{iScan} = transformStatistic(permuteFweTfceP{iScan},precision,params); 
        end
      end
      if params.bootstrapFweAdjustment
        bootstrapFweParametricP{iScan} = transformStatistic(bootstrapFweParametricP{iScan},precision,params);
        if params.TFCE 
          bootstrapFweTfceP{iScan} = transformStatistic(bootstrapFweTfceP{iScan},precision,params);
        end
      end
    end
    if params.permutationTests
      [permuteP{iScan},fdrPermuteP{iScan},fwePermuteP{iScan}] = transformStatistic(permuteP{iScan},precision,params); 
      if params.TFCE 
        [tfcePermuteP{iScan},fdrTfcePermuteP{iScan},fweTfceRandT{iScan}] = transformStatistic(tfcePermuteP{iScan},precision,params); 
      end
    end
    if params.bootstrapTests
      [bootstrapP{iScan},fdrBootstrapP{iScan},fweBootstrapP{iScan}] = transformStatistic(bootstrapP{iScan},precision,params); 
      if params.bootstrapFweAdjustment && ~params.noDoubleBootstrap
        bootstrapFweBootstrapP{iScan} = transformStatistic(bootstrapFweBootstrapP{iScan},precision,params);
      end
      if params.TFCE 
        [tfceBootstrapP{iScan},fdrTfceBootstrapP{iScan},fweTfceBootstrapP{iScan}] = transformStatistic(tfceBootstrapP{iScan},precision,params);
        if params.bootstrapFweAdjustment && ~params.noDoubleBootstrap
          bootstrapFweTfceBootstrapP{iScan} = transformStatistic(bootstrapFweTfceBootstrapP{iScan},precision,params);
        end
      end
    end
  end
  

  % save parameters
  glmAnal.d{iScan} = copyFields(d);
  stimvol = d.stimvol;
  for i=1:length(stimvol)
    stimvol{i} = unique(ceil(stimvol{i}/d.designSupersampling));
  end
  glmAnal.d{iScan}.stimvol = stimvol;
  glmAnal.d{iScan}.hrf = mrDownsample(d.hrf, d.designSupersampling/d.estimationSupersampling,floor(rem(scanParams{iScan}.acquisitionDelay,d.tr/d.estimationSupersampling)*d.designSupersampling/d.estimationSupersampling/d.tr)+1);
  glmAnal.d{iScan}.actualhrf = d.hrf;
  clear('d');
  glmAnal.d{iScan}.acquisitionDelay = scanParams{iScan}.acquisitionDelay;
  glmAnal.d{iScan}.EVnames = params.EVnames;                %this should be removed if viewGet can get params from GLM analysis
  glmAnal.d{iScan}.dim = [scanDims{iScan} numVolumes];
  glmAnal.d{iScan}.ehdr = NaN([scanDims{iScan} size(ehdr,4) size(ehdr,5)],precision);
  glmAnal.d{iScan}.ehdr(subsetBox{iScan}(1,1):subsetBox{iScan}(1,2),subsetBox{iScan}(2,1):subsetBox{iScan}(2,2),subsetBox{iScan}(3,1):subsetBox{iScan}(3,2),:,:) = ehdr;
  glmAnal.d{iScan}.s2 = NaN(scanDims{iScan},precision);
  glmAnal.d{iScan}.s2(subsetBox{iScan}(1,1):subsetBox{iScan}(1,2),subsetBox{iScan}(2,1):subsetBox{iScan}(2,2),subsetBox{iScan}(3,1):subsetBox{iScan}(3,2),:,:) = s2; %JB
  clear('ehdr','s2');
  if params.covCorrection
    glmAnal.d{iScan}.autoCorrelationParameters = NaN([scanDims{iScan} size(autoCorrelationParameters,4)],precision);
    glmAnal.d{iScan}.autoCorrelationParameters(subsetBox{iScan}(1,1):subsetBox{iScan}(1,2),subsetBox{iScan}(2,1):subsetBox{iScan}(2,2),subsetBox{iScan}(3,1):subsetBox{iScan}(3,2),:,:) = autoCorrelationParameters; %JB
    clear('autoCorrelationParameters')
  end
  if params.bootstrapTests && params.bootstrapIntervals
    glmAnal.d{iScan}.ehdrBootstrapCIs = NaN([scanDims{iScan} size(ehdrBootstrapCIs,4) size(ehdrBootstrapCIs,5)],precision);
    glmAnal.d{iScan}.ehdrBootstrapCIs(subsetBox{iScan}(1,1):subsetBox{iScan}(1,2),subsetBox{iScan}(2,1):subsetBox{iScan}(2,2),subsetBox{iScan}(3,1):subsetBox{iScan}(3,2),:,:) = ehdrBootstrapCIs; %JB
    clear('ehdrBootstrapCIs')
  end
  if numberContrasts
    glmAnal.d{iScan}.contrasts = params.contrasts;
    if params.bootstrapTests && params.bootstrapIntervals  && (nnz(params.componentsToTest)==1 || strcmp(params.componentsCombination,'Add'))
      glmAnal.d{iScan}.contrastBootstrapCIs = NaN([scanDims{iScan} size(contrastBootstrapCIs,4) size(contrastBootstrapCIs,5)],precision);
      glmAnal.d{iScan}.contrastBootstrapCIs(subsetBox{iScan}(1,1):subsetBox{iScan}(1,2),subsetBox{iScan}(2,1):subsetBox{iScan}(2,2),subsetBox{iScan}(3,1):subsetBox{iScan}(3,2),:,:) = contrastBootstrapCIs; %JB
      clear('contrastBootstrapCIs')
    end
  end
  if numberFtests
    glmAnal.d{iScan}.fTestNames = params.fTestNames;        %this should be removed if viewGet can get params from GLM analysis
    glmAnal.d{iScan}.restrictions = params.restrictions;               %this should be removed if viewGet can get params from GLM analysis
  end


end

tic
%-------------------------------------------------------- Output Analysis ---------------------------------------------------
dateString = datestr(now);
glmAnal.name = params.saveName;
if strcmp(params.hrfModel,'hrfDeconvolution')
  glmAnal.type = 'deconvAnal';
elseif isempty(params.restrictions) && ~params.computeTtests
  glmAnal.type = 'glmAnal';
else
  glmAnal.type = 'glmAnalStats';
end
glmAnal.groupName = params.groupName;
glmAnal.function = 'glmAnalysis';
glmAnal.reconcileFunction = 'defaultReconcileParams';
glmAnal.mergeFunction = 'defaultMergeParams';
glmAnal.guiFunction = 'glmAnalysisGUI';
glmAnal.params = params;
glmAnal.date = dateString;

%--------------------------------------------------------- Output overlay structures
nScans = viewGet(thisView,'nScans');
% create generic parameters 
defaultOverlay.groupName = params.groupName;
defaultOverlay.function = 'glmAnalysis';
defaultOverlay.reconcileFunction = 'defaultReconcileParams';
defaultOverlay.date = dateString;
defaultOverlay.params = cell(1,nScans);
% colormap is made with a little bit less on the dark end
defaultOverlay.colormap = hot(312);
defaultOverlay.colormap = defaultOverlay.colormap(end-255:end,:);
defaultOverlay.alpha = 1;
defaultOverlay.interrogator = 'glmPlot';
defaultOverlay.mergeFunction = 'defaultMergeParams';
defaultOverlay.colormapType = 'normal';
defaultOverlay.range = [0 1];
defaultOverlay.clip = [0 1];
defaultOverlay.alphaOverlay='';
defaultOverlay.alphaOverlayExponent=1;
defaultOverlay.data = cell(1,nScans);
defaultOverlay.name = '';
for iScan = params.scanNum
   defaultOverlay.data{iScan} = NaN(scanDims{iScan},precision); %to make values outside the box transparent
end


%------------------------------------------------------ save the r2 overlay
overlays = defaultOverlay;
overlays.name = 'r2';
overlays.colormapType = 'setRangeToMax';
for iScan = params.scanNum
   overlays.data{iScan}(subsetBox{iScan}(1,1):subsetBox{iScan}(1,2),subsetBox{iScan}(2,1):subsetBox{iScan}(2,2),subsetBox{iScan}(3,1):subsetBox{iScan}(3,2)) = r2{iScan};
   overlays.params{iScan} = scanParams{iScan};
end
clear('r2');

%--------------------------------------------- save the contrast beta weights overlay(s) (ehdr if no contrast)
contrastNames = makeContrastNames(params.contrasts,params.EVnames,params.tTestSide);
if numberContrasts && (nnz(params.componentsToTest)==1 || strcmp(params.componentsCombination,'Add'))
  %this is to mask the beta values by the probability/Z maps     
  betaAlphaOverlay = cell(numberContrasts,1);

  %find the values for the scale of the beta overlays
  thisOverlay = defaultOverlay;
  ordered_abs_betas=[];
  for iScan = params.scanNum
    ordered_abs_betas = [ordered_abs_betas; contrast{iScan}(:)];
  end
  ordered_abs_betas = ordered_abs_betas(~isnan(ordered_abs_betas));
  min_beta = min(min(min(min(min(ordered_abs_betas)))));
  max_beta = max(max(max(max(max(ordered_abs_betas)))));
  ordered_abs_betas = sort(abs(ordered_abs_betas));
  beta_perc95 = 0; 
  beta_perc95 = max(beta_perc95,ordered_abs_betas(round(numel(ordered_abs_betas)*.95))); %take the 95th percentile for the min/max
  thisOverlay.range = [-beta_perc95 beta_perc95];
  thisOverlay.clip = [min_beta max_beta];
  thisOverlay.colormap = jet(256);
  for iContrast = 1:numberContrasts
    overlays(end+1)=thisOverlay;
    overlays(end).alphaOverlay=betaAlphaOverlay{iContrast};
    overlays(end).name = contrastNames{iContrast};
    for iScan = params.scanNum
      overlays(end).data{iScan} = NaN(scanDims{iScan},precision); %to make values outside the box transparent
      overlays(end).data{iScan}(subsetBox{iScan}(1,1):subsetBox{iScan}(1,2),subsetBox{iScan}(2,1):subsetBox{iScan}(2,2),subsetBox{iScan}(3,1):subsetBox{iScan}(3,2)) =...
         contrast{iScan}(:,:,:,iContrast);
      overlays(end).params{iScan} = scanParams{iScan};
    end
  end
  clear('contrast');
end

%--------------------------------------------- save the Test overlay(s)
if numberTests
  
  testNames = [contrastNames params.fTestNames];
  for iTest = 1:numberContrasts+numberFtests
    if iTest<=numberContrasts
      testNames{iTest} = ['T (' testNames{iTest} ')'];
    else
      testNames{iTest} = ['F (' testNames{iTest} ')'];
    end
  end
  
  if params.outputStatistic
    tempStatistic = cell2mat(statistic);
    if numberContrasts && params.computeTtests
      thisOverlay = defaultOverlay;
      thisOverlay.colormap = jet(256);
      max_abs_T = max(max(max(max(max(abs(tempStatistic(:,:,:,1:numberContrasts)))))));
      thisOverlay.range = [-max_abs_T max_abs_T];
      thisOverlay.clip = thisOverlay.range;
      for iTest = 1:numberContrasts
        overlays(end+1)=thisOverlay;
        overlays(end).name = testNames{iTest};
        for iScan = params.scanNum
          overlays(end).data{iScan}(subsetBox{iScan}(1,1):subsetBox{iScan}(1,2),subsetBox{iScan}(2,1):subsetBox{iScan}(2,2),subsetBox{iScan}(3,1):subsetBox{iScan}(3,2)) =...
             statistic{iScan}(:,:,:,iTest);
          overlays(end).params{iScan} = scanParams{iScan};
        end
      end
    end
    if numberFtests
      thisOverlay = defaultOverlay;
      thisOverlay.range(1) = min(min(min(min(min(tempStatistic(:,:,:,numberContrasts+(1:numberFtests)))))));
      thisOverlay.range(2) = min(max(max(max(max(tempStatistic(:,:,:,numberContrasts+(1:numberFtests)))))));
      thisOverlay.clip = thisOverlay.range;
      for iTest = numberContrasts+(1:numberFtests)
        overlays(end+1)=thisOverlay;
        overlays(end).name = testNames{iTest};
        for iScan = params.scanNum
          overlays(end).data{iScan}(subsetBox{iScan}(1,1):subsetBox{iScan}(1,2),subsetBox{iScan}(2,1):subsetBox{iScan}(2,2),subsetBox{iScan}(3,1):subsetBox{iScan}(3,2)) =...
             statistic{iScan}(:,:,:,iTest);
          overlays(end).params{iScan} = scanParams{iScan};
        end
      end
    end
    clear('statistic','tempStatistic');
  end

  lastUncorrectedTest = 0;
  lastFDRCorrectedTest = 0;
  lastFWECorrectedTest = 0;
  if params.parametricTests
    overlays = [overlays makeOverlay(defaultOverlay, parametricP, subsetBox, params.scanNum, scanParams, ...
                                       '', params.testOutput, testNames, params.statisticalThreshold, precision)];
    lastUncorrectedTest = length(overlays);
    clear('parametricP');
    if params.bootstrapFweAdjustment
      overlays = [overlays makeOverlay(defaultOverlay, bootstrapFweParametricP, subsetBox, params.scanNum, scanParams, ...
                                         'bootstrap-FWE-adjusted ', params.testOutput, testNames, params.statisticalThreshold, precision)];
      lastFWECorrectedTest = length(overlays);
      clear('bootstrapFweParametricP');
      if params.TFCE  
        overlays = [overlays makeOverlay(defaultOverlay, bootstrapFweTfceP, subsetBox, params.scanNum, scanParams, ...
                                           'bootstrap-FWE-adjusted TFCE ', params.testOutput, testNames, params.statisticalThreshold, precision)];
        lastFWECorrectedTest = length(overlays);
        clear('bootstrapFweTfceP');
      end
    end
    if params.permutationFweAdjustment
      overlays = [overlays makeOverlay(defaultOverlay, permuteFweParametricP, subsetBox, params.scanNum, scanParams, ...
                                         'permutation-FWE-adjusted ', params.testOutput, testNames, params.statisticalThreshold, precision)];
      lastFWECorrectedTest = length(overlays);
      clear('permuteFweParametricP');
      if params.TFCE
        overlays = [overlays makeOverlay(defaultOverlay, permuteFweTfceP, subsetBox, params.scanNum, scanParams, ...
                                         'permutation-FWE-adjusted TFCE ',params.testOutput, testNames, params.statisticalThreshold, precision)];
        lastFWECorrectedTest = length(overlays);
        clear('permuteFweTfceP');
      end
    end
    if params.fdrAdjustment
      overlays = [overlays makeOverlay(defaultOverlay, fdrParametricP, subsetBox, params.scanNum, scanParams, ...
                                      'FDR-adjusted ',params.testOutput, testNames, params.statisticalThreshold, precision)];
      lastFDRCorrectedTest = length(overlays);
      clear('fdrParametricP');
    end
    if params.fweAdjustment
      overlays = [overlays makeOverlay(defaultOverlay, fweParametricP, subsetBox, params.scanNum, scanParams, ...
                                      'FWE-adjusted ',params.testOutput, testNames, params.statisticalThreshold, precision)];
      lastFWECorrectedTest = length(overlays);
      clear('fweParametricP');
    end
  end

  if params.bootstrapTests
    overlays = [overlays makeOverlay(defaultOverlay, bootstrapP, subsetBox, params.scanNum, scanParams, ...
                                       'bootstrap ', params.testOutput, testNames, params.statisticalThreshold, precision)];
    lastUncorrectedTest = length(overlays);
    clear('bootstrapP');
    if params.bootstrapFweAdjustment && ~params.noDoubleBootstrap
      overlays = [overlays makeOverlay(defaultOverlay, bootstrapFweBootstrapP, subsetBox, params.scanNum, scanParams, ...
                                         'bootstrap-FWE-adjusted bootstrap ', params.testOutput, testNames, params.statisticalThreshold, precision)];
      lastFWECorrectedTest = length(overlays);
      clear('bootstrapFweBootstrapP');
    end
    if params.fdrAdjustment
      overlays = [overlays makeOverlay(defaultOverlay, fdrBootstrapP, subsetBox, params.scanNum, scanParams, ...
                                      'FDR-adjusted bootstrap ',params.testOutput, testNames, params.statisticalThreshold, precision)];
      lastFDRCorrectedTest = length(overlays);
      clear('fdrBootstrapTp');
    end
    if params.fweAdjustment
      overlays = [overlays makeOverlay(defaultOverlay, fweBootstrapP, subsetBox, params.scanNum, scanParams, ...
                                      'FWE-adjusted bootstrap ',params.testOutput, testNames, params.statisticalThreshold, precision)];
      lastFWECorrectedTest = length(overlays);
      clear('fweBootstrapTp');
    end
    if params.TFCE  
      overlays = [overlays makeOverlay(defaultOverlay, tfceBootstrapP, subsetBox, params.scanNum, scanParams, ...
                                         'bootstrap TFCE ', params.testOutput, testNames, params.statisticalThreshold, precision)];
      lastUncorrectedTest = length(overlays);
      clear('tfceBootstrapP');
      if params.bootstrapFweAdjustment && ~params.noDoubleBootstrap
        overlays = [overlays makeOverlay(defaultOverlay, bootstrapFweTfceBootstrapP, subsetBox, params.scanNum, scanParams, ...
                                           'bootstrap-FWE-adjusted bootstrap TFCE ', params.testOutput, testNames, params.statisticalThreshold, precision)];
        lastFWECorrectedTest = length(overlays);
        clear('bootstrapFweTfceBootstrapP');
      end
      if params.fdrAdjustment
        overlays = [overlays makeOverlay(defaultOverlay, fdrTfceBootstrapP, subsetBox, params.scanNum, scanParams, ...
                                           'FDR-adjusted bootstrap TFCE ', params.testOutput, testNames, params.statisticalThreshold, precision)];
        lastFDRCorrectedTest = length(overlays);
        clear('fdrTfceBootstrapP');
      end
      if params.fweAdjustment
        overlays = [overlays makeOverlay(defaultOverlay, fweTfceBootstrapP, subsetBox, params.scanNum, scanParams, ...
                                           'FWE-adjusted bootstrap TFCE ', params.testOutput, testNames, params.statisticalThreshold, precision)];
        lastFWECorrectedTest = length(overlays);
        clear('fweTfceBootstrapP');
      end
    end
  end

  if params.permutationTests
    overlays = [overlays makeOverlay(defaultOverlay, permuteP, subsetBox, params.scanNum, scanParams, ...
                                       'permutation ', params.testOutput, testNames, params.statisticalThreshold, precision)];
    lastUncorrectedTest = length(overlays);
    clear('permuteP');
    if params.fdrAdjustment
      overlays = [overlays makeOverlay(defaultOverlay, fdrPermuteP, subsetBox, params.scanNum, scanParams, ...
                                      'FDR-adjusted permutation ',params.testOutput, testNames, params.statisticalThreshold, precision)];
      lastFDRCorrectedTest = length(overlays);
      clear('fdrPermuteP');
    end
    if params.fweAdjustment
      overlays = [overlays makeOverlay(defaultOverlay, fwePermuteP, subsetBox, params.scanNum, scanParams, ...
                                      'FWE-adjusted permutation ',params.testOutput, testNames, params.statisticalThreshold, precision)];
      lastFWECorrectedTest = length(overlays);
      clear('fwePermuteP');
    end
    if params.TFCE
      overlays = [overlays makeOverlay(defaultOverlay, tfcePermuteP, subsetBox, params.scanNum, scanParams, ...
                                      'permutation TFCE ',params.testOutput, testNames, params.statisticalThreshold, precision)];
      lastUncorrectedTest = length(overlays);
      clear('tfcePermuteP');
      if params.fdrAdjustment
        overlays = [overlays makeOverlay(defaultOverlay, fdrTfcePermuteP, subsetBox, params.scanNum, scanParams, ...
                                        'FDR-adjusted permutation TFCE ',params.testOutput, testNames, params.statisticalThreshold, precision)];
        lastFDRCorrectedTest = length(overlays);
        clear('fdrTfcePermuteP');
      end
      if params.fweAdjustment
        overlays = [overlays makeOverlay(defaultOverlay, fweTfceRandT, subsetBox, params.scanNum, scanParams, ...
                                        'FWE-adjusted permutation TFCE ',params.testOutput, testNames, params.statisticalThreshold, precision)];
        lastFWECorrectedTest = length(overlays);
        clear('fweTfceRandT');
      end
    end
  end
  
  if strcmp(params.alphaContrastOverlay,'Uncorrected') && lastUncorrectedTest
    lastContrastAlphaOverlay = lastUncorrectedTest;
  elseif strcmp(params.alphaContrastOverlay,'FDR') && lastFDRCorrectedTest
    lastContrastAlphaOverlay = lastFDRCorrectedTest;
  elseif strcmp(params.alphaContrastOverlay,'FWE') && lastFWECorrectedTest
    lastContrastAlphaOverlay = lastFWECorrectedTest;
  else
    lastContrastAlphaOverlay = 0;
  end
  if lastContrastAlphaOverlay && params.computeTtests && (nnz(params.componentsToTest)==1 || strcmp(params.componentsCombination,'Add'))
    %set the contrast alpha overlay to the statistical value
    switch(params.testOutput)
      case 'P value'                                                  %statistical maps
        betaAlphaOverlayExponent = -1;      % with inverse masking for p values
      case {'Z value','-log10(P) value'}
        betaAlphaOverlayExponent = .5;      %or normal masking for Z or log10(p) values
    end
    for iContrast = 2:numberContrasts+1
      overlays(iContrast).alphaOverlayExponent=betaAlphaOverlayExponent;
      overlays(iContrast).alphaOverlay = overlays(lastContrastAlphaOverlay-numberFtests-numberContrasts+iContrast-1).name;
    end
  end
  
end

% Try to remove unnecessary fields from the d structure in glmAnal so that it isn't too large to save
removeFields = {'stimfile',{'varname','paramInfo'}};
for iScan = 1:length(glmAnal.d)
  if ~isempty(glmAnal.d{iScan})
    for iField = 1:length(removeFields)
      % simple string, means to remove that field
      if ~iscell(removeFields{iField})
        glmAnal.d{iScan} = rmfield(glmAnal.d{iScan},removeFields{iField});
      else
        % a cell means to remove the field within a field (like d{iScan}.varname.paramInfo)
        structName = sprintf('glmAnal.d{%i}',iScan);
        for jField = 1:(length(removeFields{iField})-1)
    structName = sprintf('%s.%s',structName,removeFields{iField}{jField});
        end
        % now remove it
        eval(sprintf('%s = rmfield(%s,''%s'');',structName,structName,removeFields{iField}{end}));
      end
    end
  end
end

%-------------------------------------------------------- Set the analysis in view
glmAnal.overlays = overlays;
thisView = viewSet(thisView,'newAnalysis',glmAnal);
if ~isempty(viewGet(thisView,'fignum'))
  refreshMLRDisplay(viewGet(thisView,'viewNum'));
end
toc

%-------------------------------------------------------- Save the analysis
saveAnalysis(thisView,glmAnal.name);

oneTimeWarning('nonZeroHrfStart',0);
oneTimeWarning('tfceOutputsZeros',0);
set(viewGet(thisView,'figNum'),'Pointer','arrow');
refreshMLRDisplay(viewGet(thisView,'viewNum'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sub routines

function outputData = reshapeToRoiBox(inputData,dataPosition)
%inputData must be size ([A 1 1 B], where A equals the number of non-zero values in dataPosition
%outputData will be size ([size(dataPosition) B])
outputData = NaN([numel(dataPosition) size(inputData,4) size(inputData,5)]);
inputData = permute(inputData,[1 4 5 2 3]); %remove singleton dimensions (y and z)
outputData(dataPosition>0,:,:) = inputData;
outputData = reshape(outputData,[size(dataPosition,1) size(dataPosition,2) size(dataPosition,3) size(outputData,2) size(outputData,3)]);




function overlays = makeOverlay(overlays,data,subsetBox,scanList,scanParams,testTypeString, testOutput, testNames, threshold, precision)
  switch testOutput
    case 'P value'
      overlays.colormap = statsColorMap(256);
      namePrefix = [testTypeString 'P ['];
      overlays.range = [0 1];
      overlays.clip = [0 threshold];
    case 'Z value'
      overlays.range = [0 norminv(1-1e-16)]; %1e-16 is smallest non-zero P value output by cdf in getGlmStatistics (local functions T2p and F2p)
      overlays.clip = [norminv(1-max(threshold,1e-16)) norminv(1-1e-16)];
      namePrefix = [testTypeString 'Z ['];
    case '-log10(P) value'
      overlays.range = [0 -log10(1e-16)];
      overlays.clip = [-log10(max(threshold,1e-16)) -log10(1e-16)];
      namePrefix = ['-log10(' testTypeString 'P) ['];
  end
  overlays.range = eval([precision '(overlays.range)']);%convert range and clip to same precision as overlay 
  overlays.clip = eval([precision '(overlays.clip)']);%convert range and clip to same precision as overlay 
  overlays = repmat(overlays,1,length(testNames));
  for iTest = 1:length(testNames)
    overlays(iTest).name = [namePrefix testNames{iTest} ']'];
    for iScan = scanList
      overlays(iTest).data{iScan}(subsetBox{iScan}(1,1):subsetBox{iScan}(1,2),subsetBox{iScan}(2,1):subsetBox{iScan}(2,2),subsetBox{iScan}(3,1):subsetBox{iScan}(3,2)) =...
         data{iScan}(:,:,:,iTest);
      overlays(iTest).params{iScan} = scanParams{iScan};
    end
  end



