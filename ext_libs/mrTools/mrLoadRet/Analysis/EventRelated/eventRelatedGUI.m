% eventRelatedGUI.m
%
%      usage: eventRelatedGUI()
%         by: justin gardner
%       date: 04/05/07
%    purpose: 
%
function params = eventRelatedGUI(varargin)

% check arguments
if ~any(nargin == [0 1 2 3 4 5])
  help eventRelatedGUI
  return
end

% get the arguments
eval(evalargs(varargin));

% if called with params, then just display
if ~ieNotDefined('params')
  retval = dispParams(params);
  if isempty(retval),params = [];end
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

% set the parameter string
paramsInfo = {...
    {'groupName',groupNames,'Name of group from which to do eventRelated analysis'},...
    {'saveName','erAnal','File name to try to save as'},...
    {'applyFiltering',1,'type=checkbox','Apply the same filtering that was used before concatenation to the design matrix. This will apply the hipass filter as well as the projection analsis if they have been done by concatTSeries'},...
};

% Get parameter values
if useDefault
  params = mrParamsDefault(paramsInfo);
else
  params = mrParamsDialog(paramsInfo,'Event Related Parameters');
end

% removed this old setting since it is not going to be maintained and superceded by the standard concat
params.inplaceConcat = false;
%{'inplaceConcat',0,'type=checkbox','Concatenate all data and design matrices in memory. This runs a differrent processing stream (ask Farshad for details). If you are using a Concatenation time series do not check this.'},...

% if empty user hit cancel
if isempty(params)
  deleteView(view);
  return
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
  deleteView(view);
  return
end

% get the parameters for each scan
params.scanParams = getEventRelatedVarname(view,viewGet(view,'groupNum',params.groupName),params.scanNum,useDefault);
if isempty(params.scanParams)
  params = [];
  deleteView(view);
  return
end
% set the scan number
for i = 1:length(params.scanNum)
  params.scanParams{params.scanNum(i)}.scanNum = params.scanNum(i);
end

deleteView(view);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% just display parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = dispParams(params)

paramsInfo = {};
if isfield(params,'applyFiltering')
  paramsInfo{end+1} = {'applyFiltering',params.applyFiltering,'editable=0','Whether filtering was applied to the columns of the stimulus convolution matrix'};
end
paramsInfo{end+1} = {'analyzedScans',1,'incdec=[-1 1]',sprintf('minmax=[1 %i]',length(params.scanNum)),'editable=0'};
% get parameters for each scan
for i = 1:length(params.scanNum)
  scanNum{i} = params.scanNum(i);
  description{i} = params.scanParams{scanNum{i}}.description;
  hdrlen{i} = params.scanParams{scanNum{i}}.hdrlen;
  preprocess{i} = params.scanParams{scanNum{i}}.preprocess;
  if isfield(params.scanParams{scanNum{i}},'taskNum')
    taskNum{i} = params.scanParams{scanNum{i}}.taskNum;
  end
  if isfield(params.scanParams{scanNum{i}},'phaseNum')
    phaseNum{i} = params.scanParams{scanNum{i}}.phaseNum;
  end
  if isfield(params.scanParams{scanNum{i}},'segmentNum')
    segmentNum{i} = params.scanParams{scanNum{i}}.segmentNum;
  end
  if isfield(params.scanParams{scanNum{i}},'varname')
    varname{i} = params.scanParams{scanNum{i}}.varname;
  end
end
paramsInfo{end+1} = {'scanNum',scanNum,'group=analyzedScans','type=numeric','editable=0','scan number'};
paramsInfo{end+1} = {'tseriesFile',params.tseriesFile,'group=analyzedScans','type=string','editable=0','Name of timeseries that was analyzed'};
paramsInfo{end+1} = {'description',description,'group=analyzedScans','type=string','editable=0','Description of the analysis'};
paramsInfo{end+1} = {'hdrlen',hdrlen,'group=analyzedScans','type=numeric','editable=0','Length of response in seconds to calculate'};
paramsInfo{end+1} = {'preprocess',preprocess,'group=analyzedScans','type=string','editable=0','String of extra commands for preprocessing. Normally you will not need to set anything here, but this allows you to do corrections to the stimvols that are calculated so that you can modify the analysis. (see wiki for details)'};
if ~ieNotDefined('taskNum')
  paramsInfo{end+1} = {'taskNum',taskNum,'group=analyzedScans','type=numeric','editable=0','The task you want to use'};
end
if ~ieNotDefined('phaseNum')
  paramsInfo{end+1} = {'phaseNum',phaseNum,'group=analyzedScans','type=numeric','editable=0','The phase of the task you want to use'};
end
if ~ieNotDefined('segmentNum')
  paramsInfo{end+1} = {'segmentNum',segmentNum,'group=analyzedScans','type=numeric','editable=0','The segment of the task you want to use'};
end
if ~ieNotDefined('varname')
  paramsInfo{end+1} = {'varname',varname,'group=analyzedScans','type=string','editable=0','The variable that was analyzed'};
end
retval = mrParamsDialog(paramsInfo,'Event Related parameters');
