% getGlmEVParamsGUI.m
%
%        $Id$
%      usage: params = getGlmEVParamsGUI(thisView,params,useDefault)
%         by: julien besle
%       date: 03/12/2010
%    purpose: return EV specific parameters for GLM analysis
%

function [scanParams, params] = getGlmEVParamsGUI(thisView,params,useDefault)

keepAsking = 1;
groupNum = viewGet(thisView,'groupNum',params.groupName);
nScans = viewGet(thisView,'nScans',groupNum);
if ~isfield(params,'scanNum') || isempty(params.scanNum)
  params.scanNum = 1:nScans;
end
if isfield(params,'scanParams') && length(params.scanParams)==nScans
   scanParams = params.scanParams;
else
   % make the output as long as the number of scans
   scanParams = cell(1,nScans);
end

nStims = zeros(1,nScans);
for iScan = params.scanNum
  %get the number of events after running the pre-processing function for each scan
  disp(sprintf('Getting number of conditions from scan %d', iScan)); 
  d = loadScan(thisView, iScan, groupNum, 0);
  d = getStimvol(d,scanParams{iScan});
  d = eventRelatedPreProcess(d,scanParams{iScan}.preprocess);
  nStims(iScan) = length(d.stimvol);
  disp(sprintf('%d conditions found', nStims(iScan)));
  stimNames{iScan}=d.stimNames;
end
if ~all(nStims(params.scanNum))
  for iScan = params.scanNum
    mrWarnDlg(sprintf('No stimulus found for scan %d',iScan));
  end
  scanParams = [];
  return;
end

if ~isfield(params,'numberEVs') || isempty(params.numberEVs)
  params.numberEVs = min(nStims(params.scanNum));
end

if ~isfield(params,'numberContrasts') || isempty(params.numberContrasts)
  if isfield(params,'contrasts')
    params.numberContrasts = size(params.contrasts,1);
  elseif  strcmp(params.hrfModel,'hrfDeconvolution') 
    params.numberContrasts = 0;
  else
    params.numberContrasts = params.numberEVs;
  end
end
if ~isfield(params,'computeTtests') || isempty(params.computeTtests)
  params.computeTtests = 0;
end
if ~isfield(params,'numberFtests') || isempty(params.numberFtests)
  if isfield(params,'params') && isfield(params,'restrictions')
    params.numberFtests = length(params.restrictions);
  else
    params.numberFtests = 0;
  end
end

while keepAsking
    
  % check that the number of EVs is not greater than the smallest number of stims in a stim file
  if params.numberEVs>min(nStims(params.scanNum))
    params.numberEVs = min(nStims(params.scanNum));
    mrWarnDlg(sprintf('(getGlmEVParamsGUI) Number of EV cannot be greater the number of stimuli in the stim file (%d)',min(nStims(params.scanNum))));
  end
  %get the stimToEV matrix from each scan and make single matrix with unique stim names
  stimToEVmatrix=[];
  uniqueStimNames = [];
  for iScan = params.scanNum
    if ~isfield(scanParams{iScan}, 'stimToEVmatrix') || isempty(scanParams{iScan}.stimToEVmatrix) %|| ...
%       ~isequal(size(scanParams{iScan}.stimToEVmatrix,2),params.numberEVs)
        thisStimToEVmatrix = eye(nStims(iScan),params.numberEVs);
    else
      thisStimToEVmatrix = scanParams{iScan}.stimToEVmatrix;
    end
    %if the number of stim names has changed
    if nStims(iScan)>size(thisStimToEVmatrix,1)
      thisStimToEVmatrix=[thisStimToEVmatrix;zeros(nStims(iScan)-size(thisStimToEVmatrix,1),size(thisStimToEVmatrix,2))];
    elseif nStims(iScan)<size(thisStimToEVmatrix,1)
      thisStimToEVmatrix=thisStimToEVmatrix(1:nStims(iScan),:);
    end
    if iScan==params.scanNum(1)
      uniqueStimNames = stimNames{params.scanNum(1)};
      stimToEVmatrix = thisStimToEVmatrix;
    else
      [additionalStimNames,indices] = setdiff(stimNames{iScan},uniqueStimNames);
      uniqueStimNames = [uniqueStimNames additionalStimNames];
      stimToEVmatrix = [stimToEVmatrix;thisStimToEVmatrix(indices,:)];
    end
  end
  nUniqueStims = length(uniqueStimNames);
  
  if fieldIsNotDefined(params, 'EVnames') 
    params.EVnames = repmat({' '},1,params.numberEVs);
    if params.numberEVs <= nUniqueStims
      params.EVnames = uniqueStimNames(1:params.numberEVs);
    else
      params.EVnames(1:nUniqueStims) = uniqueStimNames;
    end
  end

  %if the number of EVs has changed
  if params.numberEVs>size(thisStimToEVmatrix,2)
    stimToEVmatrix=[thisStimToEVmatrix zeros(size(stimToEVmatrix,1),params.numberEVs-size(stimToEVmatrix,2))];
  elseif params.numberEVs<size(thisStimToEVmatrix,2)
    stimToEVmatrix=stimToEVmatrix(:,1:params.numberEVs);
  end
  if params.numberEVs>length(params.EVnames)
    for iEV=size(params.EVnames)+1:params.numberEVs
      params.EVnames{iEV}='';
    end
  elseif params.numberEVs<length(params.EVnames)
    params.EVnames = params.EVnames(1:params.numberEVs);
  end
  
  paramsInfo={};
  paramsInfo{1} = {'numberEVs',params.numberEVs,'minmax=[0 inf]','incdec=[-1 1]','incdecType=plusMinus','Number of Explanatory Variables in the model = number of columns in the design matrix. If modifying the number of F-tests, press OK to redraw the dialog with the new number of F-tests'};
  paramsInfo{2} = {'EVnames', params.EVnames, 'type=stringarray','Name of the Explanatory Variables to be estimated. EVs that are set to 0 in all rows will be removed from the model.'};
  for iEvent = 1:nUniqueStims
    paramsInfo{end+1} = {fixBadChars(uniqueStimNames{iEvent}), stimToEVmatrix(iEvent,:),...
      'incdec=[-1 1]','incdecType=plusMinus','minmax=[0 inf]',...
      ['How much stimulus ' fixBadChars(uniqueStimNames{iEvent}) ' contributes to each EV.']};
  end
        
  %get timing parameters from scan parameters
  framePeriod = viewGet(thisView,'framePeriod',params.scanNum(1),groupNum);
  supersamplingMode = 'Set value';
  designSupersampling = 1;
  estimationSupersampling = 1;
  acquisitionDelay = framePeriod/2;
  stimDuration = [];
  stimDurationMode = 'Event-related';
  durationFromFile = 1;
  for iScan = params.scanNum
    if ~fieldIsNotDefined(scanParams{iScan},'supersamplingMode')
       supersamplingMode = scanParams{iScan}.supersamplingMode;
    end
    if ~fieldIsNotDefined(scanParams{iScan},'designSupersampling')
       designSupersampling = scanParams{iScan}.designSupersampling;
    end
    if ~fieldIsNotDefined(scanParams{iScan},'estimationSupersampling')
       estimationSupersampling = scanParams{iScan}.estimationSupersampling;
    end
    if ~fieldIsNotDefined(scanParams{iScan},'acquisitionDelay')
       acquisitionDelay = scanParams{iScan}.acquisitionDelay;
    end
    if ~fieldIsNotDefined(scanParams{iScan},'stimDurationMode') 
       stimDurationMode = scanParams{iScan}.stimDurationMode;
    end
    if ~fieldIsNotDefined(scanParams{iScan},'stimDuration') 
       stimDuration = scanParams{iScan}.stimDuration;
    end
    stimfile = viewGet(thisView,'stimfile',params.scanNum(1),groupNum);
    for iFile=1:length(stimfile)
      if ~strcmp(stimfile{1}.filetype,'eventtimes') || ~isfield(stimfile{1}.mylog,'stimdurations_s')
        durationFromFile=0;
      end
    end
  end
  if durationFromFile
    stimDurationModeMenu = {'Event-related','Block','Set value','From file'};
  else
    stimDurationModeMenu = {'Event-related','Block','Set value'};
    if strcmp(stimDurationMode,'From file')
      stimDurationMode = 'Event-related';
    end
  end
  stimDurationModeMenu = putOnTopOfList(stimDurationMode,stimDurationModeMenu);

  if strcmp(params.hrfModel,'hrfDeconvolution') 
    canonicalVisibleOption = 'visible=0';
    designSupersamplingOption = 'visible=0';
    if strcmp(stimfile{1}.filetype,'eventtimes')
      deconvolutionVisibleOption = 'visible=1';
    else
      estimationSupersampling=1;
      deconvolutionVisibleOption = 'visible=0';
    end
    designSupersampling=estimationSupersampling;
    supersamplingModeMenu = {'Set value'};
  else
    estimationSupersampling=1;
    canonicalVisibleOption = 'visible=1';
    deconvolutionVisibleOption = 'visible=0';
    if strcmp(stimfile{1}.filetype,'eventtimes')
      designSupersamplingOption = 'visible=1';
    else
      designSupersamplingOption = 'visible=0';
      designSupersampling=1;
    end
    supersamplingModeMenu = {'Set value','Automatic'};
  end
  if ~strcmp(stimfile{1}.filetype,'eventtimes') %if the stimfile is not the mylog type
    stimDuration = max(framePeriod,framePeriod*round(stimDuration/framePeriod)); %round the duration to a multiple of the frame period
  end
  if ~strcmp(params.hrfModel,'hrfDeconvolution')  && strcmp(stimfile{1}.filetype,'eventtimes')
    superSamplingOption = 'visible=1';
  else
    superSamplingOption = 'visible=0';
  end
  if ismember(supersamplingMode,supersamplingModeMenu)
    supersamplingModeMenu = putOnTopOfList(supersamplingMode,supersamplingModeMenu);
  end
    
  paramsInfo = [paramsInfo {...
    {'stimDurationMode', stimDurationModeMenu, canonicalVisibleOption, 'How the duration of each stimulation is specified in the design matrix: ''Event-related'': each event has the minimum possible duration given the chosen supersampling (generally 1 TR); ''Block'': each event lasts until an event of a different type is found; ''Set value'': uses  the value specified in stimDuration parameter. ''From file'': the stimulus durations are read from the stim files (only mylog file type).'},...
    {'stimDuration', stimDuration, canonicalVisibleOption, 'type=numeric', sprintf('incdec=[-%f %f]',framePeriod,framePeriod), 'minmax=[.05 inf]','incdecType=plusMinus', 'Duration of the boxcar function that will be convolved with the HRF model at each stimulation event (in sec, resolution=0.05s). For MGL/AFNI files, the duration will be rounded to the closest multiple of the frame period (TR), excluding 0.'},...
    {'supersamplingMode', supersamplingModeMenu,superSamplingOption, 'How the design matrix super-sampling is determined. ''Set value'': supersampling value is set manually in the edit box below, usually to 1 (no super-sampling); ''Automatic'': supersampling is automatically chosen to accomodate the fine temporal resolution specified in the stim file (with a 20 Hz resolution limit).'},...
    {'estimationSupersampling',estimationSupersampling, deconvolutionVisibleOption,'incdec=[-1 1]','incdecType=plusMinus','minmax=[1 inf]','Supersampling factor of the deconvolution HRF model. Set this to more than 1 in order to resolve the estimated HDR at a temporal resolution that is less than the frame rate. This is only required if both the design and the acquisition have been designed to achieve subsample HDR estimation, and is not supported for MGL/AFNI stim files.'},...
    {'designSupersampling',designSupersampling, designSupersamplingOption,'incdec=[-1 1]','incdecType=plusMinus','minmax=[1 inf]','Supersampling factor of the GLM model. Set this to more than one in order to take into account stimulus times and duration that do not coincide with the frame period, or to round stimulus presentation times to the nearest sample onset. This is not supported for MGL/AFNI stim files.'},...
    {'acquisitionDelay',acquisitionDelay, sprintf('minmax=[0 %f]',framePeriod-0.05),'Time (in sec) at which the signal is actually acquired, on average across slices. This is normally set to half the frame period. If slice motion correction has been used, change to reference slice acquisition time. You might also need to change this value for sparse imaging designs.'},...
    {'showDesign', 0, 'type=pushbutton','buttonString=Show Design','Shows the experimental design (before convolution with an HRF model) and the design matrix.',...
        'callback',{@plotExperimentalDesign,params,scanParams,thisView,uniqueStimNames,stimNames,framePeriod},'passParams=1'}...
    {'numberContrasts',params.numberContrasts,'minmax=[0 inf]','incdec=[-1 1]','incdecType=plusMinus', 'Number of contrasts on which to perform a T-test. Both contrast values and inference test outcomes will be ouput as overlays.'},...
    {'computeTtests', params.computeTtests,'type=checkbox', 'visible=0', 'Whether contrasts are tested for significance using T-tests'},...
    {'numberFtests',params.numberFtests,'minmax=[0 inf]','incdec=[-1 1]','incdecType=plusMinus','Number of F-tests to be computed and output as overlays. An F-test tests the overall significance of a collection of contrasts. A collection of contrast is defined by a restriction matrix (one contrast per row) '},...
     }];

  % now we have all the dialog information, ask the user to set parameters
  if useDefault
     tempParams = mrParamsDefault(paramsInfo);
  else
     tempParams = mrParamsDialog(paramsInfo,'Set design parameters');
  end

  % user hit cancel
  if isempty(tempParams)
     scanParams = tempParams;
     return
  end

  currentNumberEVs = params.numberEVs;
  
  [params,scanParams]=convertStimToEVmatrix(params,tempParams,scanParams,uniqueStimNames,stimNames,framePeriod);
  
  if tempParams.numberEVs~=currentNumberEVs 
    params.numberEVs=tempParams.numberEVs;
    %if the number of EVs has changed, redraw the menu
    disp('User changed the number of EVs')
  elseif strcmp(scanParams{params.scanNum(1)}.stimDurationMode,'Set value') && isempty(scanParams{params.scanNum(1)}.stimDuration)
    mrWarnDlg('Please set a stim duration value');
  elseif length(tempParams.EVnames)~=tempParams.numberEVs
    mrWarnDlg('Please enter a name for all the EVs');
  elseif ~isequal(sort(tempParams.EVnames),unique(tempParams.EVnames))
    mrWarnDlg('All the EV names must be different');
  else
    keepAsking = 0;
  end
  
  if keepAsking && useDefault %there were incompatible parameters but this is the script mode (no GUI)
    scanParams = [];
    return;
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% convertStimToEVmatrix %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [params,scanParams]=convertStimToEVmatrix(params,tempParams,scanParams,uniqueStimNames,stimNames,framePeriod)
  % form stimToEV matrix from fields
  stimToEVmatrix = zeros(length(uniqueStimNames),params.numberEVs);
  for iEvent = 1:length(uniqueStimNames)
    stimToEVmatrix(iEvent,:) = tempParams.(fixBadChars(uniqueStimNames{iEvent}));
    tempParams = mrParamsRemoveField(tempParams,fixBadChars(uniqueStimNames{iEvent}));
  end
  EVsToRemove = find(~any(stimToEVmatrix,1));
  for iEV = EVsToRemove
    disp(sprintf('(getGlmEVParamsGUI) Removing EV ''%s'' because it is empty',tempParams.EVnames{iEV}));
  end
  stimToEVmatrix(:,EVsToRemove) = [];
  tempParams.EVnames(EVsToRemove)=[];
  %update number of EVs
  tempParams.numberEVs = size(stimToEVmatrix,2);
  
  if strcmp(params.hrfModel,'hrfDeconvolution')
    tempParams.designSupersampling = tempParams.estimationSupersampling;
  end
  if strcmp(params.hrfModel,'hrfDeconvolution') || (strcmp(tempParams.stimDurationMode,'Event-related') ...
      && ~strcmp(tempParams.supersamplingMode,'Automatic'))
    tempParams.stimDuration = framePeriod/tempParams.designSupersampling;
  elseif ismember(tempParams.stimDurationMode,{'From file','Block'})
    tempParams.stimDuration = [];
  elseif strcmp(tempParams.stimDurationMode,'Set value')
    %round the duration to a multiple of the supersampled frame period
    tempParams.stimDuration = max(framePeriod/tempParams.designSupersampling,framePeriod/tempParams.designSupersampling*round(tempParams.stimDuration/framePeriod*tempParams.designSupersampling)); 
  end

  if tempParams.numberContrasts %compute T-test for any contrast
    tempParams.computeTtests=1;
  else
    tempParams.computeTtests=0;
  end
  
  %Add stimToEVmatrix and stimNames to each scan params
  for iScan = params.scanNum
    %get subset of stimToEVmatrix for this scan
    [~,whichStims] = ismember(stimNames{iScan},uniqueStimNames);
    newParams =  mrParamsDefault({{'stimToEVmatrix',stimToEVmatrix(whichStims,:),'Matrix forming EVs from combinations of stimulus types'},...
                       {'stimNames',stimNames{iScan},'type=stringArray','Names of stimulus types'},...
                       });
    scanParams{iScan} = mrParamsCopyFields(newParams,scanParams{iScan});
    scanParams{iScan} = mrParamsCopyFields(tempParams,scanParams{iScan},{'stimDurationMode','stimDuration','supersamplingMode','designSupersampling','estimationSupersampling','acquisitionDelay'});
  end
  
  tempParams = mrParamsRemoveField(tempParams,'stimDurationMode');
  tempParams = mrParamsRemoveField(tempParams,'stimDuration');
  tempParams = mrParamsRemoveField(tempParams,'supersamplingMode');
  tempParams = mrParamsRemoveField(tempParams,'designSupersampling');
  tempParams = mrParamsRemoveField(tempParams,'estimationSupersampling');
  tempParams = mrParamsRemoveField(tempParams,'acquisitionDelay');
  tempParams = mrParamsRemoveField(tempParams,'showDesign');
  
  params = mrParamsCopyFields(tempParams,params);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% plotStimDesign %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotExperimentalDesign(thisScanParams,params,scanParams,thisView,uniqueStimNames,stimNames,framePeriod)

%we need to convert the stimToEVmatrix
[params,scanParams]=convertStimToEVmatrix(params,thisScanParams,scanParams,uniqueStimNames,stimNames,framePeriod);

if strcmp(scanParams{params.scanNum(1)}.stimDurationMode,'Set value') && isempty(scanParams{params.scanNum(1)}.stimDuration)
    mrWarnDlg('Please set a stim duration value');
    return
end
 
nScans = length(params.scanNum);
expDesignFigure = initializeFigure('plotExperimentalDesign',nScans,'horizontal');
designMatrixFigure = initializeFigure('plotDesignMatrix',nScans,'vertical');

cScan=0;
axisLength = zeros(1,length(params.scanNum));
tSeriesAxes = zeros(1,length(params.scanNum));
designMatrixAxes = zeros(1,length(params.scanNum));
for iScan = params.scanNum
  params.scanParams{iScan} = copyFields(scanParams{iScan},params.scanParams{iScan}); %we use copyFields here instead of mrParamsCopyFields because the latter only copies fields with a corresponding paramInfo entry
  colors = randomColors(length(params.EVnames));
  %replace all unused stimuli by one EV (if any)
  if any(~any(params.scanParams{iScan}.stimToEVmatrix,2))
    params.scanParams{iScan}.stimToEVmatrix(:,end+1) = ~any(params.scanParams{iScan}.stimToEVmatrix,2);
    params.EVnames{end+1} = 'Not used';
    colors(end+1,:) = [.85 .85 .85]; %last color for unused stims
  end
  
  d = loadScan(thisView, iScan, viewGet(thisView,'groupNum',params.groupName), 0);
  d = getStimvol(d,params.scanParams{iScan});
  [params.hrfParams,d.hrf] = feval(params.hrfModel, params.hrfParams, d.tr/d.designSupersampling,params.scanParams{iScan}.acquisitionDelay,1);
  d = eventRelatedPreProcess(d,params.scanParams{iScan}.preprocess);
  d = makeDesignMatrix(d,params,1,iScan);
  if  strcmp(params.EVnames{end},'Not used')
    scm = d.scm(:,1:end-d.nHrfComponents);
    EVnames = params.EVnames(1:end-1);
    nHdr = d.nhdr-1;
  else
    scm = d.scm;
    EVnames = params.EVnames;
    nHdr = d.nhdr;
  end
  testDesignMatrix(scm,nHdr,d.nHrfComponents,EVnames);
  
  cScan=cScan+1;
  
  if ~isempty(d.scm) && size(d.scm,2)==length(params.EVnames)*d.nHrfComponents
    %the length for the axis depends on the number of volumes times the frame period
    axisLength(cScan) = size(d.EVmatrix,1)/d.designSupersampling*d.tr;

    tSeriesAxes(cScan) = axes('parent',expDesignFigure,'outerposition',getSubplotPosition(1,cScan,[7 1],ones(nScans,1),0,.05));
    hold(tSeriesAxes(cScan),'on')
    [h,hTransition] = plotStims(d.EVmatrix, d.stimDurations, d.tr/d.designSupersampling, colors, tSeriesAxes(cScan),d.runTransitions);

    title(tSeriesAxes(cScan),viewGet(thisView,'description',iScan),'interpreter','none');
    %legend
    legendStrings = params.EVnames(h>0);
    h = h(h>0);
    if ~isempty(hTransition)
      h = [h hTransition];
      legendStrings = [legendStrings {'Run Transitions'}];
    end
    lhandle = legend(h,legendStrings,'position',getSubplotPosition(2,cScan,[7 1],ones(nScans,1),0,.05));
    set(lhandle,'Interpreter','none','box','off');
    
    %plot the design matrix
    designMatrixAxes(cScan) = axes('parent',designMatrixFigure,'outerposition',getSubplotPosition(cScan,1,ones(nScans,1),[7 1],0,.05));
    dimensions = size(scm);
    extended_matrix = zeros(dimensions+1);
    extended_matrix(1:size(scm,1),1:size(scm,2)) = scm;
    hMatrix = pcolor(designMatrixAxes(cScan),(1:dimensions(2)+1)-.5,(1:dimensions(1)+1)-.5,  double(extended_matrix)); %first vector must correspond to columns of the matrix
                %and second vector to the rows... (how retarded is that ?)
    set(hMatrix,'EdgeColor','none');
    Xticks = 1:size(scm,2);
    set(designMatrixAxes(cScan), 'Ydir', 'reverse');
    ylabel(designMatrixAxes(cScan),'Scans');
    %set rotated component labels
    set(designMatrixAxes(cScan), 'Xtick', Xticks );
    xTickStringsStrings = cell(1,length(EVnames)*d.nHrfComponents);
    for i = 1:length(EVnames)
      for j = 1:d.nHrfComponents
        xTickStringsStrings{(i-1)*d.nHrfComponents+j}=[EVnames{i} ' - Component ' num2str(j)];
      end
    end
    axisCoords = axis(designMatrixAxes(cScan)); % Current axis limits
    text(Xticks,axisCoords(4)*ones(1,length(Xticks)),xTickStringsStrings,...
      'parent',designMatrixAxes(cScan),'HorizontalAlignment','right',...
      'VerticalAlignment','top','Rotation',45,'interpreter','none');
    % Remove the default labels
    set(designMatrixAxes(cScan),'XTickLabel','')
    colormap(designMatrixAxes(cScan),gray);
    %add run transitions
    if ~isempty(d.runTransitions)
      hold(designMatrixAxes(cScan),'on');
      for iRun = 2:size(d.runTransitions,1)
        hTransitionDesign = plot(designMatrixAxes(cScan),axisCoords(1:2),repmat(d.runTransitions(iRun,1)/d.designSupersampling,1,2),'--r');
      end
    end
    if iScan==params.scanNum(end)
      hColorbar = colorbar('peer',designMatrixAxes(cScan));
      if size(d.runTransitions,1)>1
        colorBarPosition = get(hColorbar,'position');
        legendPosition = colorBarPosition;
        legendPosition(2)= colorBarPosition(2)/2;
        legendPosition(3)= 1-colorBarPosition(1);
        legendPosition(4)= colorBarPosition(2)/2;
        legendHandle = legend(hTransitionDesign,{'Run transitions'},'Position',legendPosition);
        set(legendHandle,'Interpreter','none','box','off');
      end
    end
    
  else
    h = axes('parent',expDesignFigure,'outerposition',getSubplotPosition(1,cScan,[7 1],ones(nScans,1),0,.05));   
    set(h,'visible','off');
    text(.5,.5,{viewGet(thisView,'description',iScan), 'Number of EVs does not match', 'Change the stimToEV matrix for this scan'},...
      'parent',h,'units','normalized','HorizontalAlignment','center');
  end
end

%resize axes according to length of each sequence
maxAxisLength = max(axisLength);
for iScan = 1:cScan
  if axisLength(iScan)
    figurePosition = get(tSeriesAxes(iScan),'position');
    figurePosition(3) = figurePosition(3)*axisLength(iScan)/maxAxisLength;
    set(tSeriesAxes(iScan),'position',figurePosition);
  end
end

function expDesignFigure = initializeFigure(figureName,nScans,mode)

expDesignFigure = selectGraphWin(0,'Make new');
set(expDesignFigure,'name',figureName);
monitorPositions = getMonitorPositions;
figurePosition = get(expDesignFigure,'position');
[whichMonitor,figurePosition]=getMonitorNumber(figurePosition,monitorPositions);
screenSize = monitorPositions(whichMonitor,:); % find which monitor the figure is displayed in
switch(mode)
  case 'horizontal'
    figurePosition(3) = screenSize(3);
    figurePosition(4) = min(screenSize(3)/8*nScans,screenSize(4));
  case 'vertical'
    figurePosition(1) = screenSize(1);
    figurePosition(2) = min(screenSize(1)/8*nScans,screenSize(2));
end
set(expDesignFigure,'position',figurePosition);
