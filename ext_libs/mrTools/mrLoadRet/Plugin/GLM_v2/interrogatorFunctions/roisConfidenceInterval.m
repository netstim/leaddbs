% roisConfidenceInterval - compute confidence interval of beta estimates in ROIs
%
%        $Id$
%      usage: [  ] = roisConfidenceInterval(thisView,overlayNum,scanNum,x,y,z,roi)
%         by: julien besle
%       date: 2010-02-15
%     inputs: 
%    outputs: 
%
%    purpose: computes confidence interval of beta estimates in ROIs
%             plots thems against numbers of events 
%        e.g:
%
function roisConfidenceInterval(thisView,overlayNum,scanNum,x,y,z,roi)

d = [];
titleString = 'RoisConfidenceInterval';

%mrGlobals;
analysisParams = viewGet(thisView,'analysisParams');
if ~isempty(analysisParams)
  analysisType = viewGet(thisView,'analysisType');
  if ismember(analysisType,{'glmAnal','glmcAnal','glmAnalStats','deconvAnal'})
    analysisParams = convertOldGlmParams(analysisParams);
    analysisParamsMenu = {'From current Analysis','Define'};
    bootstrapParams.analysisType = 'GLM';
  elseif strcmp(analysisType,'erAnal')
    analysisParamsMenu = {'From current Analysis','Define'};
    bootstrapParams.analysisType = 'Deconvolution';
  else
    analysisParams = [];
  end
end
if isempty(analysisParams)
   analysisParamsMenu = {'Define'};
   analysisParams.groupName = viewGet(thisView,'groupName',thisView.curGroup);
   bootstrapParams.analysisType = 'GLM';
end

analysisTypeMenu = {'GLM','Deconvolution'};
analysisTypeMenu = putOnTopOfList(bootstrapParams.analysisType,analysisTypeMenu);
bootstrapTypeMenu = {'Residuals','Single-trial Estimates'};   %if Deconvolution, could add third method based on modifying stimvol before call to makescm
whichROIsMenu = {'Loaded ROIs','Visible ROIs'};
bootstrapParams.maxNCombinations = Inf;
bootstrapParams.confidenceLevel = .95;
bootstrapParams.totalNRand = 1000;
maxNTimeUnits = 10; %let's not compute more than 20 different event numbers
roiAveragingMenu = {'beforeFittingGLM','afterFittingGLM'};
%analysisParams.testParams.contrasts = [];
%analysisParams.testParams.fTests = [];

askForParams = 1;
while askForParams
  paramsInfo = {...
      {'analysisType',analysisTypeMenu,'type=popupmenu','Type of analysis to estimate confidence interval for.'},...
      {'analysisParams',analysisParamsMenu,'type=popupmenu','Keep GLM parameters from saved analysis or redefine'},...
      {'bootstrapType',bootstrapTypeMenu,'type=popupmenu','At what stage bootstrap is performed'},...
      {'roiAveraging', roiAveragingMenu, 'type=popupmenu', 'averaging across voxels is performed either before or after fitting the GLM and computing the statistics'},...
      {'maxNCombinations', bootstrapParams.maxNCombinations, 'minmax=[1 inf]', 'Maximum number of combinations for a given number of runs'},...
      {'totalNRand', bootstrapParams.totalNRand, 'minmax=[1 inf]', 'Number of bootstrap resampling per combination'},...
      {'confidenceLevel', bootstrapParams.confidenceLevel, 'minmax=[0 1]', 'Defines the confidence interval boundaries. Example: if level = .95, boundaries will be at 2.5 and 97.5 % of the bootstrap distribution'},...
      {'whichROIs', whichROIsMenu, 'type=popupmenu', 'Which Rois to include in the analysis'},...
    };

  bootstrapParams = mrParamsDialog(paramsInfo, 'Bootstrap parameters');
  % Abort if params empty
  if ieNotDefined('bootstrapParams'),return,end

  analysisParamsMenu = putOnTopOfList(bootstrapParams.analysisParams,analysisParamsMenu);
  analysisTypeMenu = putOnTopOfList(bootstrapParams.analysisType,analysisTypeMenu);
  bootstrapTypeMenu = putOnTopOfList(bootstrapParams.bootstrapType,bootstrapTypeMenu);
  roiAveragingMenu = putOnTopOfList(bootstrapParams.roiAveraging,roiAveragingMenu);

  if strcmp(bootstrapParams.analysisParams,'Define')
    useDefault =0;
  else
    useDefault =1;
  end
  switch(bootstrapParams.analysisType)
     case 'GLM'
        %get GLM  parameters
        tempParams = glmAnalysisGUI('params',analysisParams,'useDefault',useDefault,'thisView',thisView);
     case 'Deconvolution'
        %get Event-related  parameters
        tempParams = eventRelatedGUI('params',analysisParams,'useDefault',useDefault);
  end
  if isempty(tempParams)
     askForParams=1;
  else
     analysisParams = tempParams;
     askForParams=0;
  end


  %controls
  switch(bootstrapParams.analysisType)
    case 'GLM'
       % add the constraint that contrasts and fTests is not compatible with Single-trial Estimates bootstrapping 
       if strcmp(bootstrapParams.bootstrapType,'Single-trial Estimates') && (~isempty(analysisParams.testParams.fTests) || ~isempty(analysisParams.testParams.contrasts))
          mrWarnDlg('''Single-trial Estimates'' not compatible with contrasts or fTests, please redefine parameters');
          askForParams=1;
       end
       if ischar(analysisParams.scanParams{params.scanNum(1)}.stimDuration) && strcmp(analysisParams.scanParams{params.scanNum(1)}.stimDuration,'fromFile')
          mrWarnDlg('Function not implemented for stimDuration = fromFile, Please redefine parameters');
          askForParams=1;
       end
       if analysisParams.covCorrection
          mrWarnDlg('Function not implemented for Covariance Correction, Please redefine parameters');
          askForParams=1;
       end

    case 'Deconvolution'
  %          % add the constraint that Single-trial Estimates bootstrapping  is not compatible with deconvolution
  %          if strcmp(bootstrapParams.bootstrapType,'Single-trial Estimates') && (~isempty(analysisParams.testParams.fTests) || ~isempty(analysisParams.testParams.contrasts))
  %             mrWarnDlg('''Single-trial Estimates'' not compatible with deconvolution analysis, please redefine parameters');
  %             askForParams=1;
  %          end
       if isfield(analysisParams.scanParams{1},'inplaceConcat') && analysisParams.scanParams{1}.inplaceConcat
          mrWarnDlg('''inplaceConcat'' is not supported, Please redefine parameters');
          askForParams=1;
       end
  end
end


alpha = (1-bootstrapParams.confidenceLevel);
maxNCombinations = bootstrapParams.maxNCombinations;
totalNRand = bootstrapParams.totalNRand;
if strcmp(bootstrapParams.analysisType,'GLM') && strcmp(analysisParams.testParams.tTestSide,'Both')
   t_alpha = alpha/2;
else
   t_alpha = alpha;
end

% set the group
thisView = viewSet(thisView,'groupName',analysisParams.groupName);

switch(bootstrapParams.whichROIs)
  case 'Loaded ROIs'
    roiNums = 1:length(thisView.ROIs);
  case 'Visible ROIs'
    roiNums = viewGet(thisView,'visibleROIs');
end
ROIs = thisView.ROIs(roiNums);
nRois = length(roiNums);
roiLabels=cell(nRois,1);
for iRoi = 1:nRois
  roiLabels{iRoi} = {ROIs(iRoi).name, ['(' num2str(size(ROIs(iRoi).coords,2)) ' voxels)']};
end


%--------------------------------------------------------- Main loop over scans ---------------------------------------------------
set(viewGet(thisView,'figNum'),'Pointer','watch');drawnow;
for scanNum = analysisParams.scanNum
   titleString = [titleString ' - ' viewGet(thisView,'description')];

   [d, roiVoxelIndices] = loadScanRois(thisView,scanNum,roiNums);
   
   switch(bootstrapParams.roiAveraging)
      case 'beforeFittingGLM'
         %average over ROIs
         actual_data = NaN(nRois,1,1,d.dim(4));
         for iRoi = 1:nRois
            actual_data(iRoi,:) = mean(d.data(roiVoxelIndices{iRoi},:,:,:),1);
         end
         d.data = actual_data; 
         d.dim = size(d.data);
   end
   
   nVoxels = size(d.data,1);
   d.volumes = 1:d.dim(4);
   
   % get the stim volumes, if empty then abort
%    if strcmp(bootstrapParams.analysisType,'GLM')
%       d.stimDuration = analysisParams.scanParams{scanNum}.stimDuration;
%       d.forceStimOnTR = analysisParams.scanParams{scanNum}.forceStimOnTR;
%    end
   d = getStimvol(d,analysisParams.scanParams{scanNum});
   if isempty(d.stimvol),mrWarnDlg('No stim volumes found');return,end
   
   switch(bootstrapParams.analysisType)
   case 'GLM'
      titleString = [titleString ' - GLM'];      
      fTests = analysisParams.testParams.fTests;
      contrasts = analysisParams.testParams.contrasts;
      %-------------- get the HRF model --------------
      [analysisParams.hrfParams,d.hrf] = feval(analysisParams.hrfModel, analysisParams.hrfParams, d.tr/d.designSupersampling,0,1);
      nEstimates = size(d.hrf,2);
      
   case 'Deconvolution'
      titleString = [titleString ' - Deconvolution'];      
      fTests = [];
      contrasts = [];
      nEstimates = ceil(analysisParams.scanParams{scanNum}.hdrlen/d.tr);
   end
   
   
   % do any call for preprocessing
   if ~isempty(analysisParams.scanParams{scanNum}.preprocess)
      d = eventRelatedPreProcess(d,analysisParams.scanParams{scanNum}.preprocess);
   end

   if isempty(contrasts) && isempty(fTests)
      nEstimateTypes = length(d.stimvol);
      plotThreshold = 0;
      estimateTypes = 'Beta Estimates';
   else
      estimateTypes = 'Contrasts/F-tests';
      nEstimateTypes = size(contrasts,1)+size(fTests,1);
      nEstimates = 1;                 
      plotThreshold = 1;
   end

   
%--------------------------------------------------------- Bootstrap ---------------------------------------------------------

   h_wait = waitbar(0,['Bootstrap Resampling for scan ' num2str(scanNum)]);
   
   switch(bootstrapParams.bootstrapType)
      
      case 'Residuals'    %BOOTSTRAP BEFORE FITTING MODEL (ON RESIDUALS)
         %compute the design matrix for all events
         runTransition = d.concatInfo.runTransition;
         switch(bootstrapParams.analysisType)
            case 'GLM'
               d = makeDesignMatrix(d,analysisParams,1,1);
            case 'Deconvolution'
               d = makescm(d,nEstimates,analysisParams.applyFiltering);
         end
         runNumbers = 1:size(runTransition,1);
         timeUnit = 'runs';
         titleString = [titleString ' - Residuals Bootstrapping'];
      
      case 'Single-trial Estimates'  %BOOTSTRAP AFTER FITTING SINGLE-EVENT MODEL 
         
         titleString = [titleString ' - Single Response Estimates Bootstrapping'];
         %first create design matrix with one single event per column
         boot_d = d;
         boot_d.stimvol = num2cell(cell2mat(d.stimvol));
         boot_d.nhdr = length(d.stimvol);
         %compute the design matrix
         switch(bootstrapParams.analysisType)
            case 'GLM'
               boot_d = makeDesignMatrix(boot_d,analysisParams,1,1);
            case 'Deconvolution'
               boot_d = makescm(boot_d,nEstimates,analysisParams.applyFiltering);
         end
         %fit the single-event model
         waitbar(0,h_wait,['Fitting Single-Event Mode l for scan ' num2str(scanNum)]);
         boot_d = getGlmStatistics(boot_d);

         %number of event per type
         for i_event =  1:nEstimateTypes  
            nEvent(i_event) = length(d.stimvol{i_event});
         end
                  %plot difference between single-event GLM Model and whole GLM model, see end of file

         %equalize the number of event per type and average across rois
         nEventPerType = min(nEvent);  %number of events per type
         min_nEventPerType = 2; %bootstrp doesn't seem to deal very well with having only one element on the first dimension of this_ehdr
         runStep = ceil((nEventPerType-min_nEventPerType)/maxNTimeUnits);   
         runNumbers = min_nEventPerType + (rem((nEventPerType-min_nEventPerType),runStep):runStep:nEventPerType-min_nEventPerType);
         timeUnit = 'events';
        
         nEstimateTypes = length(nEvent);
         ehdr = NaN(nEventPerType,nEstimateTypes,nEstimates,nRois);
         for iRoi = 1:nRois
            c_type = 0;
            for i_type =  1:nEstimateTypes  %for each event-type
               switch(bootstrapParams.roiAveraging)
                  case 'beforeFittingGLM'
                     ehdr(:,i_type,:,iRoi) = boot_d.ehdr(iRoi,1,1,c_type+1:c_type+nEventPerType,:);
                  case 'afterFittingGLM'
                     ehdr(:,i_type,:,iRoi) = squeeze(mean(boot_d.ehdr(roiVoxelIndices{iRoi},1,1,c_type+1:c_type+nEventPerType,:),1));
               end
               c_type = c_type + nEvent(i_type);
            end
         end
   end
   
   actualValue = zeros(nEstimateTypes, nEstimates, nRois,length(runNumbers));
   bootstrapMean = actualValue;
   lowerErrorBar = actualValue;
   upperErrorBar = actualValue;
   if isempty(contrasts) && isempty(fTests)
    statisticalThreshold = [];
   else
    statisticalThreshold = actualValue;
   end
   
   maxNCombinations = min(maxNCombinations,totalNRand);
   waitbar(0,h_wait,['Bootstrap Resampling for scan ' num2str(scanNum)]);
   accumulated_runs = 0;
 
   for iRun = 1:length(runNumbers)   
      runLabels{iRun} = [num2str(runNumbers(iRun)) ' ' timeUnit];
      fprintf(1, ['(roisConfidenceInterval) ' runLabels{iRun} '...']);
      
      if runNumbers(end)<20 %if the number of runs is small enough,
         % get all the possible combinations from the current number of runs
         combinations =  combnk(1:runNumbers(end),runNumbers(iRun));
         if size(combinations,1) > maxNCombinations %and take 
            %take a random subset of all the combinations
            permuted_indices = randperm(size(combinations,1));
            permuted_indices = permuted_indices(1:maxNCombinations);
            combinations = combinations(permuted_indices,:);
         end
      else %if number of runs is >20, take combinations at random. it is highly unlikely to take twice the same anyway
         combinations = NaN(min(maxNCombinations,totalNRand),runNumbers(end));
         for i_combi = 1:min(maxNCombinations,totalNRand)
            combinations(i_combi,:) = randperm(runNumbers(end));
         end
         combinations = combinations(:,1:runNumbers(iRun));
      end
      
      n_combinations(iRun) = size(combinations,1);
      nBootstrap(iRun) = round(totalNRand/n_combinations(iRun));
      actual_n_bootstrap(iRun) = nBootstrap(iRun)*n_combinations(iRun);
      lower_index = max(1,round((alpha/2)*nBootstrap(iRun)));
      upper_index = round((1-alpha/2)*nBootstrap(iRun));
      
      tic
      for i_combi = 1:n_combinations(iRun)
         waitbar((accumulated_runs + ((i_combi-1)/n_combinations(iRun))*iRun)  / sum(cumsum(ones(1,length(runNumbers)))),h_wait);
         
         switch(bootstrapParams.bootstrapType)

            case 'Single-trial Estimates'  %AFTER-FITTING SINGLE-EVENT MODEL STUFF
               this_ehdr = ehdr(1:runNumbers(iRun),:,:,:);
               actualValue(:,:,:,iRun) = actualValue(:,:,:, iRun) + shiftdim(mean(this_ehdr,1),1);
               bootstrap_ehdr = bootstrp(nBootstrap(iRun),@(x)mean(x,1),this_ehdr);
               bootstrap_ehdr = reshape(bootstrap_ehdr,nBootstrap(iRun),nEstimateTypes,nEstimates,nRois);
              
               bootstrap_ehdr = sort(bootstrap_ehdr,1);
               bootstrapMean(:,:,:,iRun) = bootstrapMean(:,:,:,iRun) + shiftdim(mean(bootstrap_ehdr,1),1);
               upperErrorBar(:,:,:,iRun) = upperErrorBar(:,:,:,iRun) + shiftdim(bootstrap_ehdr(upper_index,:,:,:) - mean(bootstrap_ehdr,1),1);
               lowerErrorBar(:,:,:,iRun) = lowerErrorBar(:,:,:,iRun) + shiftdim(mean(bootstrap_ehdr,1) - bootstrap_ehdr(lower_index,:,:,:),1);

            case 'Residuals'  %BEFORE-FITTING MODEL STUFF
               %take the corresponding subset of data and stimvol
               boot_d = d;
               samples = [];
               for j_run = 1:size(combinations,2)
                  samples = [samples runTransition(combinations(i_combi,j_run),1):runTransition(combinations(i_combi,j_run),2)];
               end
               boot_d.data = d.data(:,:,:,samples);
               boot_d.scm = d.scm(samples,:);
               boot_d.nFrames = length(samples);
               boot_d.volumes = 1:boot_d.nFrames;
               boot_d.dim(4) = boot_d.nFrames;

               % Compute parameter estimates
               analysisParams.nBootstrap = nBootstrap(iRun);
               [boot_d, dump, T, F] = getGlmStatistics(boot_d, analysisParams,0);
               if ~isempty(T)
                  T = reshape(T(:,1,1,:,:),nVoxels,size(contrasts,1),nEstimates,nBootstrap(iRun));
               end
               if ~isempty(F)
                  F = reshape(F(:,1,1,:,:),nVoxels,size(fTests,1),nEstimates,nBootstrap(iRun));
               end

               %remove y and z dimensions (=1) 
               if isempty(contrasts) && isempty(fTests)
                  bootstrap_values = reshape(boot_d.ehdr(:,1,1,:,:,:),nVoxels,nEstimateTypes,nEstimates,nBootstrap(iRun));
               else
                  bootstrap_values = cat(2,T,F);
                  statisticalThreshold(1:size(contrasts,1),:,:,iRun) = icdf('T',1-t_alpha,boot_d.rdf);
                  for i_f = 1: size(fTests,1)
                     statisticalThreshold(size(contrasts,1) + i_f,:,:,iRun) = icdf('f', 1-alpha, boot_d.mdf(i_f), boot_d.rdf);
                  end
               end
               %put the voxel dimension at the end
               bootstrap_values = permute(bootstrap_values,[2 3 4 1]);
               clear('boot_d');

               %we will estimate the extent of the CI per combination, and then get an average over combinations
               for iRoi = 1:nRois
                  %find bootstrap mean and bottom and upper bound of the CI
                  switch(bootstrapParams.roiAveraging)
                     case 'beforeFittingGLM'
                        ehdr_roi = bootstrap_values(:,:,:,iRoi);
                     case 'afterFittingGLM'
                        %here we want an estimation for the CI *for the roi*, so we need to compute an average of the value for each roi *before* sorting
                        ehdr_roi = mean(bootstrap_values(:,:,:,roiVoxelIndices{iRoi}),4);
                  end
                  actualValue(:,:,iRoi,iRun) = actualValue(:,:,iRoi,iRun)+ehdr_roi(:,:,1);
                  ehdr_roi = sort(ehdr_roi,3);
                  roi_mean = mean(ehdr_roi,3);
                  upperErrorBar(:,:,iRoi,iRun) = upperErrorBar(:,:,iRoi,iRun) + ehdr_roi(:,:,upper_index) - roi_mean;
                  lowerErrorBar(:,:,iRoi,iRun) = lowerErrorBar(:,:,iRoi,iRun) + roi_mean - ehdr_roi(:,:,lower_index);
                  bootstrapMean(:,:,iRoi,iRun) = bootstrapMean(:,:,iRoi,iRun) + roi_mean;
               end
         end
         
      end
      toc
      actualValue(:,:,:,iRun) = actualValue(:,:,:,iRun)/n_combinations(iRun);
      upperErrorBar(:,:,:,iRun) = upperErrorBar(:,:,:,iRun)/n_combinations(iRun);
      lowerErrorBar(:,:,:,iRun) = lowerErrorBar(:,:,:,iRun)/n_combinations(iRun);
      bootstrapMean(:,:,:,iRun) = bootstrapMean(:,:,:,iRun)/n_combinations(iRun);
      
      accumulated_runs = accumulated_runs + iRun;
   end
   n_bootstrap_string = {['Number of combination of ' timeUnit ':' num2str(n_combinations)],...
                         ['Number of resamples per combination: ' num2str(nBootstrap)],...
                         ['Total number of resamples: ' num2str(actual_n_bootstrap)]};
                      
   close(h_wait);
 
   
   %---------------------- construct figure for this scan -----------------%
   
   selectGraphWin;
   mrGlobals
   figNum = MLR.graphFigure;
   set(figNum,'Name',titleString);

   for i_est = 1:nEstimateTypes
      if isempty(contrasts) || i_est<=size(contrasts,1) 
         if isempty(contrasts) 
            event_number = i_est;
         else
            event_number = find(contrasts(i_est,:)~=0);
         end
         if length(event_number)==1
            if event_number<=length(d.stimNames)
               contrast_name = d.stimNames{event_number};
            else
               contrast_name = ['event type ' num2str(event_number)];
            end
         else
             contrast_name = num2str(contrasts(i_est,:));
         end
         if isempty(contrasts)
            estimateTypeNames{i_est} = [num2str(i_est) ': ' contrast_name];
         else
            estimateTypeNames{i_est} = [num2str(i_est) ': T (' contrast_name ')'];
         end
         estimateTypeNum(i_est) = i_est;
      else
         estimateTypeNames{i_est} = [num2str(i_est+1) ': F (' num2str(fTests(i_est-size(contrasts,1),:)) ')'];
         estimateTypeNum(i_est) = i_est+1;
      end
   end
   
  %put ROIs dimension last, makes things easier after
  actualValue = permute(actualValue,[1 2 4 3]);
  bootstrapMean = permute(bootstrapMean,[1 2 4 3]);
  lowerErrorBar = permute(lowerErrorBar,[1 2 4 3]);
  upperErrorBar = permute(upperErrorBar,[1 2 4 3]);
  statisticalThreshold = permute(statisticalThreshold,[1 2 4 3]);
  
  switch(bootstrapParams.analysisType)
    case 'GLM'
      plotMode = 'Estimate Histograms';
      estimateComponentsUnit = 'Components';
      estimateComponentsName = 'Basis Functions';
      estimateComponents = 1:nEstimates;

    case 'Deconvolution'
      plotMode = 'Estimate Curves';
      estimateComponentsUnit = 'Time (sec)';
      estimateComponentsName = 'TRs';
      estimateComponents = ((1:nEstimates)-.5)*d.tr;
  end
  componentPrefix = 'Component ';
  componentLabels = mat2cell(sprintf([componentPrefix '%2.d'],estimateComponents),1,(length(componentPrefix)+2)*ones(1,length(estimateComponents)));

  axesLabels = {estimateTypeNames componentLabels runLabels roiLabels};
  axesUnits = {estimateTypes estimateComponentsUnit timeUnit 'ROIs'};
  axesValues = {estimateTypeNum estimateComponents runNumbers 1:nRois};
  axesNames = {estimateTypes, estimateComponentsName, 'Runs','ROIs'};
  
  uicontrol('Parent',figNum, 'String',n_bootstrap_string,'style','text','unit','normalized','position',[.1 .01  .8 .06]);
  plotConfidenceIntervals(actualValue,bootstrapMean,lowerErrorBar,upperErrorBar,statisticalThreshold,...
    axesNames,axesValues,axesUnits,axesLabels,figNum,plotMode);
  
end

end



%%%%%%%%----------- Old
% to plot difference between single-event GLM Model and whole GLM model
%                         %fit the original model
%                         waitbar(0,h_wait,['Fitting Original Model for scan ' num2str(scanNum)]);
%                         d = getGlmStatistics(d);
%
%                         %compute parameter estimates over ROIs
%                         for iRoi = 1:nRois
%                            d_ehdr_roi(:,iRoi) = squeeze(mean(d.ehdr(roiVoxelIndices{iRoi},1,1,:)));
%                            c_event = 0;
%                            for i_event =  1:length(d.stimvol)  %for each event-type
%                               boot_d_ehdr_roi(i_event,iRoi) = mean(squeeze(mean(boot_d.ehdr(roiVoxelIndices{iRoi},1,1,c_event+1:c_event+nEvent(i_event)))));
%                               c_event = c_event + nEvent(i_event);
%                            end
%                         end
% 
%                         %compare parameter estimates
%                         figure;
%                         figureData.minValue = min(min(d_ehdr_roi(:)),min(boot_d_ehdr_roi(:)));
%                         figureData.maxValue = max(max(d_ehdr_roi(:)),max(boot_d_ehdr_roi(:)));
%                         for iRoi = 1:nRois
%                            roi_names{iRoi} = ROIs(iRoi).name;
%                         end
% 
%                         subplot(1,2,1)
%                         bar(d_ehdr_roi');
%                         set(gca,'Xticklabel',roi_names);
%                         scale = axis;
%                         scale([3 4]) = [figureData.minValue figureData.maxValue];
%                         axis(scale);
%                         legend('EV1', 'EV2', 'EV3', 'EV4', 'EV5' );
%                         title(['Original Model. Mean RSS = ' num2str(mean(d.rss)) ' - Mean R2 = ' num2str(mean(d.r2))]);
% 
%                         subplot(1,2,2)
%                         bar(boot_d_ehdr_roi');
%                         set(gca,'Xticklabel',roi_names);
%                         scale = axis;
%                         scale([3 4]) = [figureData.minValue figureData.maxValue];
%                         axis(scale);
%                         legend('EV1', 'EV2', 'EV3', 'EV4', 'EV5' );
%                         title(['Single-event Model. Mean RSS = ' num2str(mean(boot_d.rss)) ' - Mean R2 = ' num2str(mean(boot_d.r2))]);





