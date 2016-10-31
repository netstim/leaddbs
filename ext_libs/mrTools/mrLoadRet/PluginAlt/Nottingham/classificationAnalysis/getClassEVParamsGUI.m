% getClassEVParamsGUI.m
%
%        $Id: getClassEVParamsGUI.m 2699 2013-03-18 16:11:42Z julien $
%      usage: params = getClassEVParamsGUI(thisView,params,useDefault)
%         by: julien besle
%       date: 03/12/2010
%    purpose: return EV specific parameters for ER/searchlight classification analysis, saved from getGlmEVParamsGUI.m (27/03/2013)
%               because original function has became incompatible with cassification analysis
%

function [scanParams, params] = getClassEVParamsGUI(thisView,params,useDefault)

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
if ~isfield(params, 'EVnames') || ~isequal(length(params.EVnames),params.numberEVs) 
  EVnames = {};
else
  EVnames = params.EVnames;
end

while keepAsking
  nStims = zeros(1,nScans);
  for iScan = params.scanNum
    %get the number of events after running the pre-processing function for each scan
    disp(sprintf('Getting number of conditions from scan %d', iScan)); 
    d = loadScan(thisView, iScan, groupNum, 0);
    d = getStimvol(d,scanParams{iScan});
    d = eventRelatedPreProcess(d,scanParams{iScan}.preprocess);
    nStims(iScan) = length(d.stimvol);
    disp(sprintf('%d conditions found', nStims));
    stimNames{iScan}=d.stimNames;
  end
  
  if ~params.numberEVs
    params.numberEVs = max(nStims);
  end
    
  %get the stimToEV matrix from each scan and make single matrix with unique stim names
  stimToEVmatrix=[];
  uniqueStimNames = [];
  for iScan = params.scanNum
    if ~isfield(scanParams{iScan}, 'stimToEVmatrix') || isempty(scanParams{iScan}.stimToEVmatrix) || ...
      ~isequal(size(scanParams{iScan}.stimToEVmatrix,2),params.numberEVs)
        thisStimToEVmatrix = eye(nStims(iScan),params.numberEVs);
    else
      thisStimToEVmatrix = scanParams{iScan}.stimToEVmatrix;
    end
    if nStims(iScan)>size(thisStimToEVmatrix,1)
      thisStimToEVmatrix=[thisStimToEVmatrix;zeros(nStims(iScan)-size(thisStimToEVmatrix,1),size(thisStimToEVmatrix,2))];
    elseif nStims(iScan)<size(thisStimToEVmatrix,1)
      thisStimToEVmatrix=thisStimToEVmatrix(1:nStims(iScan),:);
    end
    if iScan==params.scanNum(1)
      uniqueStimNames = stimNames{1};
      stimToEVmatrix = thisStimToEVmatrix;
    else
      [additionalStimNames,indices] = setdiff(stimNames{iScan},uniqueStimNames);
      uniqueStimNames = [uniqueStimNames additionalStimNames];
      stimToEVmatrix = [stimToEVmatrix;thisStimToEVmatrix(indices,:)];
    end
  end
  nUniqueStims = length(uniqueStimNames);
  
  if isempty(EVnames) 
    EVnames = repmat({' '},1,params.numberEVs);
    if params.numberEVs <= nUniqueStims
      EVnames = uniqueStimNames(1:params.numberEVs);
    else
      EVnames(1:nUniqueStims) = uniqueStimNames;
    end
  end

  paramsInfo{1}= {'EVnames', EVnames, 'type=stringarray','Name of the Explanatory Variables to be estimated. EVs that are set to 0 in all rows will be removed from the model.'};
  for iEvent = 1:nUniqueStims
    paramsInfo{iEvent+1} = {fixBadChars(uniqueStimNames{iEvent}), stimToEVmatrix(iEvent,:),...
      'incdec=[-1 1]','incdecType=plusMinus','minmax=[0 inf]',...
      ['How much stimulus ' fixBadChars(uniqueStimNames{iEvent}) ' contributes to each EV.']};
  end
  paramsInfo{end+1}= {'showDesign', 0, 'type=pushbutton','buttonString=Show Design','Shows the experimental design (before convolution with an HRF model) and the design matrix.)',...
          'callback',{@plotExperimentalDesign,params,scanParams,thisView,uniqueStimNames,stimNames},'passParams=1'};

  %%%%%%%%%%%%%%%%%%%%%%%
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

  [params,scanParams]=convertStimToEVmatrix(params,tempParams,scanParams,uniqueStimNames,stimNames);
  
  if keepAsking && useDefault %there were incompatible parameters but this is the script mode (no GUI)
    scanParams = [];
    return;
  end
%   end
  keepAsking = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% convertStimToEVmatrix %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [params,scanParams]=convertStimToEVmatrix(params,tempParams,scanParams,uniqueStimNames,stimNames)
  % form stimToEV matrix from fields
  stimToEVmatrix = zeros(length(uniqueStimNames),params.numberEVs);
  for iEvent = 1:length(uniqueStimNames)
    stimToEVmatrix(iEvent,:) = tempParams.(fixBadChars(uniqueStimNames{iEvent}));
  end
  EVnames = tempParams.EVnames;
  EVsToRemove = find(~any(stimToEVmatrix,1));
  for iEV = EVsToRemove
    disp(sprintf('(getGlmEVParamsGUI) Removing EV ''%s'' because it is empty',EVnames{iEV}));
  end
  stimToEVmatrix(:,EVsToRemove) = [];
  EVnames(EVsToRemove)=[];
  %update number of EVs
  params.numberEVs = size(stimToEVmatrix,2);
  
  %Add stimToEVmatrix and stimNames to each scan params
  for iScan = params.scanNum
    %get subset of stimToEVmatrix fro this scan
    [~,whichStims] = ismember(stimNames{iScan},uniqueStimNames);
    newParams =  mrParamsDefault({{'stimToEVmatrix',stimToEVmatrix(whichStims,:),'Matrix forming EVs from combinations of stimulus types'},...
                       {'stimNames',stimNames{iScan},'type=strinArray','Names of stimulus types'}});
    scanParams{iScan} = mrParamsCopyFields(newParams,scanParams{iScan});
  end
  params.EVnames = EVnames;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% plotStimDesign %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotExperimentalDesign(thisScanParams,params,scanParams,thisView,uniqueStimNames,stimNames)

%we need to convert the stimToEVmatrix
[params,scanParams]=convertStimToEVmatrix(params,thisScanParams,scanParams,uniqueStimNames,stimNames);
 
nScans = length(scanParams);
expDesignFigure = initializeFigure('plotExperimentalDesign',nScans,'horizontal');
designMatrixFigure = initializeFigure('plotDesignMatrix',nScans,'vertical');

cScan=0;
axisLength = zeros(1,length(params.scanNum));
tSeriesAxes = zeros(1,length(params.scanNum));
designMatrixAxes = zeros(1,length(params.scanNum));
for iScan = params.scanNum
  params.scanParams{iScan} = copyFields(scanParams{iScan},params.scanParams{iScan}); %we use copyFields here instead of mrParamsCopyFields because the latter only copies fields with a corresponding paramInfo entry
  EVnames = thisScanParams.EVnames;
  colors = randomColors(length(EVnames));
  %replace all unused stimuli by one EV (if any)
  if any(~any(params.scanParams{iScan}.stimToEVmatrix,2))
    params.scanParams{iScan}.stimToEVmatrix(:,end+1) = ~any(params.scanParams{iScan}.stimToEVmatrix,2);
    EVnames{end+1} = 'Not used';
    colors(end+1,:) = [.85 .85 .85]; %last color for unused stims
  end
  
  d = loadScan(thisView, iScan, viewGet(thisView,'groupNum',params.groupName), 0);
  d = getStimvol(d,params.scanParams{iScan});
  [params.hrfParams,d.hrf] = feval(params.hrfModel, params.hrfParams, d.tr/d.designSupersampling,params.scanParams{iScan}.acquisitionDelay,1);
  d = eventRelatedPreProcess(d,params.scanParams{iScan}.preprocess);
  d = makeDesignMatrix(d,params,1,iScan);
  cScan=cScan+1;
  
  if ~isempty(d.scm) && size(d.scm,2)==length(EVnames)*d.nHrfComponents
    %the length for the axis depends on the number of volumes times the frame period
    axisLength(cScan) = size(d.EVmatrix,1)/d.designSupersampling*d.tr;

    tSeriesAxes(cScan) = axes('parent',expDesignFigure,'outerposition',getSubplotPosition(1,cScan,[7 1],ones(nScans,1),0,.05));
    hold(tSeriesAxes(cScan),'on')
    [h,hTransition] = plotStims(d.EVmatrix, d.stimDurations, d.tr/d.designSupersampling, colors, tSeriesAxes(cScan),d.runTransitions);

    title(tSeriesAxes(cScan),viewGet(thisView,'description',iScan),'interpreter','none');
    %legend
    legendStrings = EVnames(h>0);
    h = h(h>0);
    if ~isempty(hTransition)
      h = [h hTransition];
      legendStrings = [legendStrings {'Run Transitions'}];
    end
    lhandle = legend(h,legendStrings,'position',getSubplotPosition(2,cScan,[7 1],ones(nScans,1),0,.05));
    set(lhandle,'Interpreter','none','box','off');
    
    %plot the design matrix
    designMatrixAxes(cScan) = axes('parent',designMatrixFigure,'outerposition',getSubplotPosition(cScan,1,ones(nScans,1),[7 1],0,.05));
    dimensions = size(d.scm);
    extended_matrix = zeros(dimensions+1);
    extended_matrix(1:size(d.scm,1),1:size(d.scm,2)) = d.scm;
    hMatrix = pcolor(designMatrixAxes(cScan),(1:dimensions(2)+1)-.5,(1:dimensions(1)+1)-.5,  double(extended_matrix)); %first vector must correspond to columns of the matrix
                %and second vector to the rows... (how retarded is that ?)
    set(hMatrix,'EdgeColor','none');
    Xticks = 1:size(d.scm,2);
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
