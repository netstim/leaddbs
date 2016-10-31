function view = loadAnalysis(view,pathStr,startPathStr,name)
%
% view = loadAnalysis(view,[pathStr],[startPathStr],[name])
%
% Loads an analysis and adds it to view.analyses.
%
% If pathStr is not specified, prompts user to select a file. If pathStr is
% specified, it is assumed to be relative to startPathStr which is assumed
% (if not specified) to be the current group's data directory. pathStr can
% be a string specifying a single analysis file or it can be a cell array
% of filenames to load multiple analyses at once.
%
% name: optional replacement name(s) for the analysis when loaded. If
% pathStr is a cell array, then name must be a cell array of the same
% length.
%
% The file must contain a structure or structures, each with the following
% fields:
% - name: string
% - groupName: string group name
% - function: string function name by which it was computed
% - reconcileFunction: string function name that reconciles params
%   and data with tseries file names.
%      [newparams,newdata] = reconcileFunction(groupName,params,data)
% - guiFunction: string function name that allows user to specify
%   or change params.
%      params = guiFunction('groupName',groupName,'params',params)
% - params: structure specifying arguments to function
%      To recompute: view = function(view,params)
% - overlays: struct array of overlays
% - curOverlay: integer corresponding to currently selected overlay
% - date: specifies when it was computed
%
% In addition the 'name' field is set to the variable name to ensure
% consistency.
%
% djh, 1/9/98
% 6/2005, djh, update to mrLoadRet-4.0

% check arguments
if ~any(nargin == [1:4])
  help loadAnalysis
  return
end

mrGlobals

% Path to analyses
if ieNotDefined('startPathStr')
  startPathStr = viewGet(view,'dataDir');
end

% Complete pathStr
if ieNotDefined('pathStr')
  pathStr = mlrGetPathStrDialog(startPathStr,'Choose one or more analyses','*.mat','on');
else
  if iscell(pathStr)
    for p=1:length(pathStr)
      pathStr{p} = fullfile(startPathStr,pathStr{p});
      pathStr{p} = [stripext(pathStr{p}),'.mat'];
    end
  else
    pathStr = {[stripext(fullfile(startPathStr,pathStr)),'.mat']};
  end
end
if isempty(pathStr)
  return
end

if ieNotDefined('name')
  name = [];
elseif ~iscell(name)
  name = {name};
end

if ~iscell(pathStr)
  pathStr = {pathStr};
end

% Load the file. Loop through the variables that were loaded and add
% each of them as a new analysis, setting analysis.fieldnames as we go.
for p = 1:length(pathStr)
  if ~exist(pathStr{p},'file')
    oldPathStr = pathStr{p};
    % if file doesn't exist, try filename/filename (e.g. erAnal is usually in erAnal/erAnal.mat)
    pathStr{p} = fullfile(fileparts(pathStr{p}),stripext(getLastDir(pathStr{p})),getLastDir(pathStr{p}));
    % if still doesn't exist, revert
    if ~exist(pathStr{p},'file'),pathStr{p} = oldPathStr;end
  end
  if exist(pathStr{p},'file')
    h = mrMsgBox(['Loading analysis: ',pathStr{p},'. Please wait']);
    s = load(pathStr{p});
    varNames = fieldnames(s);
    analysis = eval(['s.',varNames{1}]);
    if isempty(name)
      analysis.name = varNames{1};
    else
      analysis.name = name{p};
    end
    % Add it to the view
    view = viewSet(view,'newAnalysis',analysis);
    mrCloseDlg(h);
  else
    mrWarnDlg(['Analysis ',pathStr{p},' not found.']);
  end
end

return;

% Test/debug
view = loadAnalysis(MLR.views{1},'co');
view = loadAnalysis(MLR.views{1},{'co','amp','ph'});
view = loadAnalysis(MLR.views{1});

