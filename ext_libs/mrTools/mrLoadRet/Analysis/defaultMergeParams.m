% defaultMergeParams.m
%
%      usage: defaultMergeParams(groupName,oldParams,newParams,oldData,newData)
%         by: justin gardner
%       date: 05/21/07
%    purpose: default function merge parameters. called by
%    saveAnalysis and saveOverlay
%
function [mergedParams mergedData] = defaultMergeParams(groupName,oldParams,newParams,oldData,newData)

% check arguments
if ~any(nargin == [3 5])
  help defaultMergeParams
  return
end

% default arguments
if ieNotDefined('oldData')
    oldData = [];
end
if ieNotDefined('newData')
    newData = [];
end

mergedParams = newParams;
mergedData = newData;

% deal with cell array of parameters
if iscell(newParams)
  for i = 1:length(oldParams)
    if isempty(newParams{i}) && ~isempty(oldParams{i})
      mergedParams{i} = oldParams{i};
      if length(oldData) >= i
	mergedData{i} = oldData{i};
      end
    end
  end
  return
end

% get scan numbers
if isfield(oldParams,'scanList') && isfield(newParams,'scanList')
  scanListName = 'scanList';
elseif isfield(oldParams,'scanNum') && isfield(newParams,'scanNum')
  scanListName = 'scanNum';
else
  scanListName = '';
end

% get list of fields to copy
scanFields = {};
if isfield(oldParams,'scanParams') && isfield(newParams,'scanParams')
  scanFields{end+1} = 'scanParams';
end

if ~isempty(scanListName)
  % go through the scan list of the old params
  for i = 1:length(oldParams.(scanListName))
    % get this scan number
    thisScanNum = oldParams.(scanListName)(i);
    % check to see if it exist in the newParams
    if isempty(find(mergedParams.(scanListName) == thisScanNum))
      % add the scan number to the merged params if it doesn't
      mergedParams.(scanListName)(end+1) = thisScanNum;
      % and tseies filename
      mergedParams.tseriesFile{end+1} = oldParams.tseriesFile{i};
      % add the fields that need to be copied
      for iFields = 1:length(scanFields)
	mergedParams.(scanFields{iFields})(thisScanNum) = ...
	    oldParams.(scanFields{iFields})(thisScanNum);
      end
      % and copy the data
      if ~isempty(newData)
	if length(oldData)>=thisScanNum
	  mergedData{thisScanNum} = oldData{thisScanNum};
	else
	  mergedData{thisScanNum} = [];
	end
      end
    end
  end
end

