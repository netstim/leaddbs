% concatTSeriesReconcileParams.m
%
%      usage: concatTSeriesReconcileParams()
%         by: justin gardner
%       date: 10/12/06
%    purpose: 
%
function params = concatTSeriesReconcileParams(groupName,params)

% check arguments
if ~any(nargin == [1 2])
  help concatTSeriesReconcileParams
  return
end


% generate params
if ieNotDefined('params')
  params.groupName = groupName;
  groupNum = viewGet([],'groupNum',groupName);
  if ~isfield(params,'scanList')
    view = newView;
    view = viewSet(view, 'groupName', groupName)
    params.scanList = selectInList(view,'scans');
  end
  params.newGroupName = 'Concatenation';
  params.baseScan = params.scanList(1);
  params.description = sprintf('Concatenation from %s of scans: %s',params.groupName,num2str(params.scanList));
  % set to one to use highpass filter
  params.filterType = 'Detrend and highpass';
  params.filterCutoff = 0.01;
  % set to one to convert to percent signal change
  params.percentSignal = 1;
  params.warp = 0;
  params.interpMethod = 'nearest';
end

