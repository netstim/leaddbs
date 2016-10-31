% eventRelatedGlmGUI.m
%
%      usage: eventRelatedGlmGUI()
%         by: farshad moradi
%       date: 06/14/07
%    purpose: 
%
function params = eventRelatedGlmGUI(varargin)

% check arguments
if ~any(nargin == [0 1 2 3 4 5])
  help eventRelatedGUI
  return
end

% get the arguments
eval(evalargs(varargin));

% if called with params, then just display
if ~ieNotDefined('params')
  dispParams(params);
  return
end

% get a view
view = newView;

% get the group names
if ieNotDefined('groupName')
  groupNames = viewGet(view,'groupNames');
else
  % if passed in name, put that on top of list to make it the default
  groupNames = putOnTopOfList(groupName,viewGet(view,'groupNames'));
end

% check for variable to just useDefaults rather than bring up gui
if ieNotDefined('useDefault')
  useDefault = 0;
end


% check if some parameters are alreay set
if ieNotDefined('params')
  params.saveName = 'glmAnal';
  params.hrfModel ='hrfDiffGamma';
  params.trSupersampling = 1;
end

if ~isfield(params, 'contrast')
  params.contrast = {'None','All','User defined'};
end

askForParams = 1;

while askForParams
    % set the parameter string
    paramsInfo = {...
        {'groupName',groupNames,'Name of group from which to do eventRelated analysis'},...
        {'saveName',params.saveName,'File name to try to save as'},...
        {'hrfModel',params.hrfModel,'Name of the function that defines the hrf used in glm'},...
        {'trSupersampling', params.trSupersampling, 'minmax=[1 100]', 'TR supersampling factor (1=no supersampling) reulting design matrix will be downsampled afterwards'},...
        {'contrast', params.contrast, 'Set to none to compute r2 overlay only. Enter all to compute an overlay for each condition (betaWeight map). Set to User defined if you want to compute a custom contrast - like the difference between condition 1 and condition 4 for example - this will bring up another dialog box where you will set the contrast of interest.'},...
    };

    % Get parameter values
    if useDefault
      params = mrParamsDefault(paramsInfo);
    else
      params = mrParamsDialog(paramsInfo, 'GLM parameters');
    end

    % if empty user hit cancel
    if isempty(params),return,end

    if strcmp(params.contrast,'User defined');
      contrastParams = mrParamsDialog({{'contrast',[],'Vector defining the contrast of interest. For instance to get a contrast between condtion 1 and 3 (when you have four conditions) you would enter [1 0 -1 0]'}},'Input contrast to compute');
      if isempty(contrastParams),return,end
      params.contrast = contrastParams.contrast;
    end
    
    % get hrf model parameters
    hrfParamsInfo = feval(params.hrfModel, 'params');

    if useDefault
      params.hrfParams = mrParamsDefault(hrfParamsInfo);
    else
      params.hrfParams = mrParamsDialog(hrfParamsInfo, sprintf('hrfModel:%s', params.hrfModel));
    end
    % if empty user hit cancel, go back
    if isempty(params.hrfParams) && ~useDefault
        askForParams = 1;
    else
        askForParams = 0;
    end
end


if ~isfield(params.hrfParams, 'description')
    params.hrfParams.description = params.hrfModel;
end

% get scans
view = viewSet(view,'groupName',params.groupName);
if ~ieNotDefined('scanList')
  params.scanNum = scanList;
elseif useDefault
  params.scanNum = 1:viewGet(view,'nScans');
else
  params.scanNum = selectInList(view,'scans');
end
if isempty(params.scanNum)
  params = [];
  return
end

% get the parameters for each scan
params.scanParams = getEventRelatedParams(view,params,useDefault);
if isempty(params.scanParams)
  params = [];
  return
end
% set the scan number
for i = 1:length(params.scanNum)
  params.scanParams{params.scanNum(i)}.scanNum = params.scanNum(i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to get the variable name that the user wants
% to do the event related analysis on, puts up a gui
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function scanParams = getEventRelatedParams(view,params,useDefault)

% make the output as long as the number of scans
scanParams = cell(1,viewGet(view,'nScans',viewGet(view,'groupNum',params.groupName)));

% check for stimfile, and if it is mgl/type then ask the
% user which variable they want to do the anlysis on
for scanNum = 1:length(params.scanNum)
  % get scan and default description
  scanInfo = sprintf('%i: %s',params.scanNum(scanNum),viewGet(view,'description',params.scanNum(scanNum)));
  description = sprintf('Event related analysis of %s: %i',params.groupName,params.scanNum(scanNum));
  % standard parameters to set
  taskVarParams = {...
      {'scan',scanInfo,'type=statictext','Description of scan to set parameters for (not editable)'},...
      {'description',description,'Event related analysis of [x...x]','Description of the analysis'}...
      {'hdrlen',25,'Length of response in seconds to calculate'}...
      {'preprocess','','String of extra commands for preprocessing. Normally you will not need to set anything here, but this allows you to do corrections to the stimvols that are calculated so that you can modify the analysis. (see wiki for details)'}...
		  };

  % make sure we are running on a set with a stimfile
  stimfile = viewGet(view,'stimfile',params.scanNum(scanNum));
  
  if isempty(stimfile)
    mrMsgBox(sprintf('No associated stimfile with scan %i in group %s',params.scanNum(scanNum),params.groupName));
    scanParams = [];
    return
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % see if we have a stimfile from mgl, in which case we should
  % ask the user what the variable name is that they want ot use for the analysis
  if strfind(stimfile{1}.filetype,'mgl')

    % check to see what style this is, if the task variable does
    % not have a segmentTrace then it mus be an old style, in which
    % we used channels
    task = cellArray(stimfile{1}.task,2);
    if isfield(stimfile{1}.myscreen,'traces') && ~isfield(task{1}{1},'segmentTrace')
      % this is the old style, get the stimtrace number
      taskVarParams{end+1} = {'stimtrace',stimfile{1}.myscreen.stimtrace,'the trace number that contains the stimulus','incdec=[-1 1]',sprintf('minmax=[%i %i]',stimfile{1}.myscreen.stimtrace,size(stimfile{1}.myscreen.traces,1))};
    else
      if exist('getTaskVarnames') ~= 2
	mrErrorDlg('(eventRelatedGUI) MGL function getTaskVarnames is not in path. You must have mgl in the path to extract stimulus timing information from an mgl stim file');
      end
      % this is the new tyle, ask for a variable name
      [varnames varnamesStr] = getTaskVarnames(stimfile{1}.task);
      % if there is more than one task, then ask the user for that
      task = cellArray(stimfile{1}.task,2);
      if length(task)>1
	taskVarParams{end+1} = {'taskNum',num2cell(1:length(task)),'The task you want to use'};
      end
      % if there are multiple phases, then ask for that
      maxPhaseNum = 0;
      maxSegNum = 0;
      for tnum = 1:length(task)
	phaseNum{tnum} = num2cell(1:length(task{tnum}));
	maxPhaseNum = max(maxPhaseNum,length(task{tnum}));
	% if there are multiple _segments_, then ask for that
	for pnum = 1:length(task{tnum})
	  segNum{tnum}{pnum} = num2cell(1:length(task{tnum}{pnum}.segmin));
	  maxSegNum = max(maxSegNum,length(segNum{tnum}{pnum}));
	end
      end
      if maxPhaseNum > 1
	if length(task) == 1
	  taskVarParams{end+1} = {'phaseNum',phaseNum{1},'The phase of the task you want to use'};
	else
	  taskVarParams{end+1} = {'phaseNum',phaseNum,'The phase of the task you want to use','contingent=taskNum'};
	end
      end
      
      % if there is more than one segement in any of the phases, ask the user to specify
      % should add some error checking.
      taskVarParams{end+1} = {'segmentBegin',1,'The segment of the trial from which your stimulus started','incdec=[-1 1]',sprintf('minmax=[0 %i]',maxSegNum)};
      taskVarParams{end+1} = {'segmentEnd',maxSegNum,'The segment of the trial at which your stimulus ended','incdec=[-1 1]',sprintf('minmax=[0 %i]',maxSegNum)};
      
      % set up to get the variable name from the user
      taskVarParams{end+1} ={'varname',varnames{1},sprintf('Analysis variables: %s',varnamesStr)};
    end
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  % give the option to use the same variable for all
  if (scanNum == 1) && (length(params.scanNum)>1)
    taskVarParams{end+1} = {'sameForAll',1,'type=checkbox','Use the same variable name for all analyses'};
  end
  %%%%%%%%%%%%%%%%%%%%%%%
  % now we have all the dialog information, ask the user to set parameters
  if useDefault
    scanParams{params.scanNum(scanNum)} = mrParamsDefault(taskVarParams);
  else
    scanParams{params.scanNum(scanNum)} = mrParamsDialog(taskVarParams);
  end
  % user hit cancel
  if isempty(scanParams{params.scanNum(scanNum)})
    scanParams = [];
    return
  end
  %%%%%%%%%%%%%%%%%%%%%%%
    
  % check if the varname is a cell array, then convert to a cell array
  % instead of a string this is so that the user can specify a variable
  % name like {{'varname'}}
  if (isfield(scanParams{params.scanNum(scanNum)},'varname') &&...
      isstr(scanParams{params.scanNum(scanNum)}.varname) && ...
      (length(scanParams{params.scanNum(scanNum)}.varname) > 1) && ...
      (scanParams{params.scanNum(scanNum)}.varname(1) == '{'))
    scanParams{params.scanNum(scanNum)}.varname = eval(scanParams{params.scanNum(scanNum)}.varname);
  end

  % if sameForAll is set, copy all parameters into all scans and break out of loop
  if isfield(scanParams{params.scanNum(scanNum)},'sameForAll') && ...
	scanParams{params.scanNum(scanNum)}.sameForAll
    for i = 2:length(params.scanNum)
      % set the other scans params to the same as this one
      scanParams{params.scanNum(i)} = scanParams{params.scanNum(1)};
      % change the description field appropriately for this scan num
      description = scanParams{params.scanNum(1)}.description;
      groupNameLoc = strfind(description,params.groupName);
      if ~isempty(groupNameLoc)
	description = sprintf('%s%s: %i',description(1:groupNameLoc(1)),params.groupName,params.scanNum(i));
      end
      scanParams{params.scanNum(i)}.description = description;
    end
    break
  end
  taskVarParams = {};
end

