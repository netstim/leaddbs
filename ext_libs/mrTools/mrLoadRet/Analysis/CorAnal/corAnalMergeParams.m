function [mergedParams,mergedData] = corAnalMergeParams(groupName,...
  oldParams,newParams,oldData,newData)
%
% [mergedParams,mergedData] = corAnalMergeParams(groupName,oldParams,newParams,[oldData],[newData])
%
% Merges corAnal params. Called by saveOverlay or saveAnalysis when file
% already exists and 'merge' option is selected. Calls
% corAnalReconcileParams to make sure that the merged params and data are
% consistent with the tseries files.
%
% oldParams: corAnal params structure in saved file.
% newParmas: corAnal params structure corresponding to newly computed
%            analysis.
% oldData: cell array of data arrays corresponding to oldParams.
%            Default: []
% newData: cell array of data arrays corresponding to newParams.
%            Default: []
%
% mergedParams: uses newParams if recompute param is non-zero. Otherwise,
%            uses oldParams.
% mergedData: corresponding merged data.
%
% djh 5/2007

if ieNotDefined('oldData')
  oldData = [];
end
if ieNotDefined('newData')
  newData = [];
end

% Reconcile params and data with current tseries and put the old and new
% params and data in the same order corresponding to that of the current
% tseries/scans. [This is now done by saveOverlay and saveAnalysis]
%[oldParams,oldData] = corAnalReconcileParams(groupName,oldParams,oldData);
%[mergedParams,mergedData] = corAnalReconcileParams(groupName,newParams,newData);
mergedParams = newParams;
mergedData = newData;

% Get nScans
groupNum = viewGet([],'groupNum',groupName);
nScans = viewGet([],'nscans',groupNum);

for scan = 1:nScans
  if ~newParams.recompute(scan)
    mergedParams.recompute(scan) = oldParams.recompute(scan);
    mergedParams.ncycles(scan) = oldParams.ncycles(scan);
    mergedParams.detrend(scan) = oldParams.detrend(scan);
    mergedParams.spatialnorm(scan) = oldParams.spatialnorm(scan);
    mergedParams.tseriesfile(scan) = oldParams.tseriesfile(scan);
    if iscell(oldData)
      mergedData{scan} = oldData{scan};
    end
  end
end
