% eventRelatedSeriesReconcileParams.m
%
%      usage: eventRelatedReconcileParams()
%         by: justin gardner
%       date: 10/20/06
%    purpose: 
%
function [params,data] = eventRelatedReconcileParams(groupName,params,data)

% check arguments
if ~any(nargin == [1 2 3 4])
  help eventRelatedReconcileParams
  return
end

% generate params if they do not exist
if ieNotDefined('params') || ~isfield(params,'groupName')
  params.groupName = groupName;
end
params.groupNum = viewGet([],'groupNum',groupName);
if ~isfield(params,'scanNum')
  view = newView;
  params.scanNum = selectInList(view,'scans');
end
if ~isfield(params,'description')
  params.description = sprintf('Event related analysis of scan %i in group %s',params.scanNum,params.groupName);
end
if ~isfield(params,'hdrlen')
  params.hdrlen = 25;
end
if ~isfield(params,'preprocess')
  params.preprocess = '';
end

