% defaultMergeParams.m
%
%      usage: defaultMergeParams(groupName,oldParams,newParams,oldData,newData)
%         by: justin gardner
%       date: 05/21/07
%    purpose: default function merge parameters. called by
%    saveAnalysis and saveOverlay
%
function [mergedParams mergedData] = pRFMergeParams(groupName,oldParams,newParams,oldData,newData)

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
  % merge any overalys - allowing nan points to be
  % overwritten by whoever has data (this allows an
  % old overlay in which partial information was
  % calculated to be merged into new overlay data
  for iData = 1:length(mergedData)
    % see if we have both old and new data
    if (length(newData) >= iData) && (length(oldData) >= iData) && ~isempty(newData{iData}) && ~isempty(oldData{iData}) && isnumeric(newData{iData}) && isnumeric(oldData{iData})
      % oldData has points not in merged data
      oldDataNotInMerged = find(~isnan(oldData{iData}) & isnan(mergedData{iData}));
      mergedData{iData}(oldDataNotInMerged) = oldData{iData}(oldDataNotInMerged);
      % newData has points not in merged data
      newDataNotInMerged = find(~isnan(newData{iData}) & isnan(mergedData{iData}));
      mergedData{iData}(newDataNotInMerged) = newData{iData}(newDataNotInMerged);
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
    else
      % does exist, so merge the two, get which points are missing
      [dump missingPoints] = setdiff(oldData{thisScanNum}.linearCoords,newData{thisScanNum}.linearCoords);
      % grab old and new linear coords and make sure that they are both row
      % vectors
      oldLinearCoords = oldData{thisScanNum}.linearCoords(missingPoints);
      oldLinearCoords = oldLinearCoords(:)';
      newLinearCoords = newData{thisScanNum}.linearCoords;
      newLinearCoords = newLinearCoords(:)';
      % and combine
      mergedData{thisScanNum}.linearCoords = [newLinearCoords oldLinearCoords];
      mergedData{thisScanNum}.params = [newData{thisScanNum}.params oldData{thisScanNum}.params(:,missingPoints)];
      mergedData{thisScanNum}.r = [newData{thisScanNum}.r' oldData{thisScanNum}.r(missingPoints,:)']';
    end
  end
end

