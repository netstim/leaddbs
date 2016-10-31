% getGlmEVParamsGUI.m
%
%        $Id: getGlmEVParamsGUI.m 2040 2011-02-13 12:20:51Z julien $
%      usage: params = getGlmEVParamsGUI(thisView,params,useDefault)
%         by: julien besle
%       date: 03/12/2010
%    purpose: return EV specific parameters for GLM analysis
%

function [scanParams, params] = getClassEventParamsGUI(thisView,params,useDefault)

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
if ~isfield(params, 'EVnames') || ~isequal(length(params.EVnames),params.numberEvents) 
  EVnames = {};
else
  EVnames = params.EVnames;
end

while keepAsking
  % check for stimfile, and if it is mgl/type then ask the
  % user which variable they want to do the anlysis on
  for iScan = params.scanNum
    
    %get the number of events after running the pre-processing function
    disp(sprintf('Getting number of conditions from scan %d', iScan)); 
    d = loadScan(thisView, iScan, groupNum, 0);
    d = getStimvol(d,scanParams{iScan});
    d = eventRelatedPreProcess(d,scanParams{iScan}.preprocess);
    nStims = length(d.stimvol);
    disp(sprintf('%d conditions found', nStims));
    
    if ~params.numberEvents
      params.numberEvents = nStims;
    end
    if ~isfield(scanParams{iScan}, 'stimToEVmatrix') || isempty(scanParams{iScan}.stimToEVmatrix) || ...
      ~isequal(size(scanParams{iScan}.stimToEVmatrix,1),nStims) || ~isequal(size(scanParams{iScan}.stimToEVmatrix,2),params.numberEvents)
        scanParams{iScan}.stimToEVmatrix = eye(nStims,params.numberEvents);
    end
    if isempty(EVnames) 
      EVnames = repmat({' '},1,params.numberEvents);
      if params.numberEvents < nStims
        EVnames = d.stimNames(1:params.numberEvents);
      else
        EVnames(1:params.numberEvents) = d.stimNames;
      end
    end

    paramsInfo = cell(1,nStims+3);
    paramsInfo{1}= {'scan', scanParams{iScan}.scan, 'editable=0','description of the scan'};
    paramsInfo{2}= {'EVnames', EVnames, 'type=stringarray','Name of the Events to be classified'};
    for iEvent = 1:nStims
      paramsInfo{iEvent+2} = {fixBadChars(d.stimNames{iEvent}), scanParams{iScan}.stimToEVmatrix(iEvent,:),...
        'incdec=[-1 1]','incdecType=plusMinus','minmax=[0 1]',...
        ['How much stimulus ' fixBadChars(d.stimNames{iEvent}) ' contributes to each Event']};
    end
%     paramsInfo{nStims+3}= {'showDesign', 0, 'type=pushbutton','buttonString=Show Experimental Design','Shows the experimental design (before convolution iwht the HRF model)',...
%             'callback',{@plotExperimentalDesign,scanParams,params,iScan,thisView,d.stimNames},'passParams=1'};

    % give the option to use the same variable for remaining scans
    if (iScan ~= params.scanNum(end)) && (length(params.scanNum)>1)
      paramsInfo{nStims+3} = {'sameForNextScans',iScan == params.scanNum(1),'type=checkbox','Use the same parameters for all scans'};
    else
      paramsInfo(nStims+3) = [];
    end

    %%%%%%%%%%%%%%%%%%%%%%%
    % now we have all the dialog information, ask the user to set parameters
    if useDefault
       tempParams = mrParamsDefault(paramsInfo);
    else
       tempParams = mrParamsDialog(paramsInfo,'Set Stimulus/EV matrix');
    end

    % user hit cancel
    if isempty(tempParams)
       scanParams = tempParams;
       return
    end
    
%     tempParams = mrParamsRemoveField(tempParams,'showDesign');
    
    % form stimToEV matrix from fields
    stimToEVmatrix = zeros(nStims,params.numberEvents);
    for iEvent = 1:nStims
      stimToEVmatrix(iEvent,:) = tempParams.(fixBadChars(d.stimNames{iEvent}));
      tempParams = mrParamsRemoveField(tempParams,fixBadChars(d.stimNames{iEvent}));
    end
    scanParams{iScan} = mrParamsCopyFields(...
      mrParamsDefault({{'stimToEVmatrix',stimToEVmatrix,'Matrix forming EVs from combinations of stimulus types'}}),scanParams{iScan});
    EVnames = tempParams.EVnames;
%     tempParams = mrParamsRemoveField(tempParams,'EVnames');
    
    scanParams{iScan} = mrParamsCopyFields(tempParams,scanParams{iScan});

    % if sameForNextScans is set, copy all parameters into all remaining scans and break out of loop
    if isfield(scanParams{iScan},'sameForNextScans') && ...
       scanParams{iScan}.sameForNextScans
       for jScan = params.scanNum(find(params.scanNum>iScan,1,'first'):end)
          % set the other scans params to the same as this one
          scanParams{jScan} = mrParamsCopyFields(scanParams{jScan},scanParams{iScan});
       end
       break
    end
  %          paramsInfo = {};
    if keepAsking && useDefault %there were incompatible parameters but this is the script mode (no GUI)
      scanParams = [];
      return;
    end
  end
  params.EVnames = EVnames;
  keepAsking = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% plotStimDesign %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotExperimentalDesign(thisScanParams,scanParams,params,scanNum,thisView,stimNames)

%we need to convert the stimToEVmatrix
for iEvent = 1:length(stimNames)
  thisScanParams.stimToEVmatrix(iEvent,:) = thisScanParams.(fixBadChars(stimNames{iEvent}));
end

if ~isfield(thisScanParams,'sameForNextScans') || thisScanParams.sameForNextScans
  %if it's the last scan or sameForNextScans is selected, we show all scans
  %copy params to remaining scans
  for jScan = params.scanNum(find(params.scanNum==scanNum,1,'first'):end)
    scanParams{jScan} = copyFields(thisScanParams,scanParams{jScan});
  end
  scanList = params.scanNum;
else
  %otherwise we only show the current scan
  scanParams = cell(1,scanNum);
  scanParams{scanNum} = thisScanParams;
  scanList = scanNum;
end

fignum = selectGraphWin(0,'Make new');
set(fignum,'name','plotExperimentalDesign');
nScans = length(scanList);
monitorPositions = getMonitorPositions;
figurePosition = get(fignum,'position');
[whichMonitor,figurePosition]=getMonitorNumber(figurePosition,monitorPositions);
screenSize = monitorPositions(whichMonitor,:); % find which monitor the figure is displayed in
figurePosition(3) = screenSize(3);
figurePosition(4) = min(screenSize(3)/8*nScans,screenSize(4));
set(fignum,'position',figurePosition);

cScan=0;
axisLength = zeros(1,length(params.scanNum));
tSeriesAxes = zeros(1,length(params.scanNum));
thisScanParams.EVnames{end+1} = 'Not used';
for iScan = scanList
  params.scanParams{iScan} = copyFields(scanParams{iScan},params.scanParams{iScan}); %we use copyFields here instead of mrParamsCopyFields because the latter only copies fields with a corresponding paramInfo entry
  %replace all unused stimuli by one EV
  params.scanParams{iScan}.stimToEVmatrix(:,end+1) = ~any(params.scanParams{iScan}.stimToEVmatrix,2);
  
  d = loadScan(thisView, iScan, viewGet(thisView,'groupNum',params.groupName), 0);
  d = getStimvol(d,params.scanParams{iScan});
  [params.hrfParams,d.hrf] = feval(params.hrfModel, params.hrfParams, d.tr/d.designSupersampling,0,1);
  d = eventRelatedPreProcess(d,params.scanParams{iScan}.preprocess);
  d = makeDesignMatrix(d,params,1,iScan);
  cScan=cScan+1;
  
  if ~isempty(d.scm) && size(d.scm,2)==length(thisScanParams.EVnames)*d.nHrfComponents
    %the length for the axis depends on the number of volumes times the frame period
    axisLength(cScan) = size(d.EVmatrix,1)*d.tr;
    colors = randomColors(length(thisScanParams.EVnames));
    colors(end,:) = [.85 .85 .85]; %last color for unused stims

    tSeriesAxes(cScan) = axes('parent',fignum,'outerposition',getSubplotPosition(1,cScan,[7 1],ones(nScans,1),0,.05));
    hold on
    [h,hTransition] = plotStims(d.EVmatrix, d.stimDurations, d.tr/d.designSupersampling, colors, tSeriesAxes(cScan),d.runTransitions);

    title(viewGet(thisView,'description',iScan),'interpreter','none');
    %legend
    legendStrings = thisScanParams.EVnames(h>0);
    h = h(h>0);
    if ~isempty(hTransition)
      h = [h hTransition];
      legendStrings = [legendStrings {'Run Transitions'}];
    end
    lhandle = legend(h,legendStrings,'position',getSubplotPosition(2,cScan,[7 1],ones(nScans,1),0,.05));
    set(lhandle,'Interpreter','none','box','off');
  else
    h = axes('parent',fignum,'outerposition',getSubplotPosition(1,cScan,[7 1],ones(nScans,1),0,.05));   
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
