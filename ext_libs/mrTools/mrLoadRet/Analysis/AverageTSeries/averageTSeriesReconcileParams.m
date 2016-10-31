function newparams = averageTSeriesReconcileParams(groupName,params,varargin)
% params = averageTSeriesReconcileParams(groupName,[params])
%
% Checks for consistency between averageTSeries parameters and current
% status of group (nscans and tseries filenames). Reconciles params by
% matching tseries filenames.
%
% params: Optional initial parameters (see averageTSeries).
% default if not passed: initialize new params with default values.
%
% djh 7/2006

groupNum = viewGet([],'groupNum',groupName);
nScans = viewGet([],'nscans',groupNum);

% evaluate other arguments
eval(evalargs(varargin));

% Get tseries filenames for this group
tseriesfiles = cell(1,nScans);
for scan = 1:nScans
    tseriesfiles{scan} = viewGet([],'tseriesFile',scan,groupNum);
end

if ieNotDefined('params')
    % Use default params
    if ieNotDefined('scanList')
      scanList = [1:nScans];
    end
    newparams.scanList = scanList;
    newparams.shiftList = zeros(size(scanList));
    newparams.reverseList = zeros(size(scanList));
    newparams.tseriesfiles = {};
    for i = 1:length(scanList)
      newparams.tseriesfiles{end+1} = tseriesfiles{scanList(i)};
    end
    newparams.baseScan = 1;
    newparams.groupName = groupName;
    newparams.aveGroupName = 'Averages';
    newparams.interpMethod = 'nearest';
    newparams.description = ['Average from ',groupName,' of scans: ',num2str(scanList)];
    newparams.fileName = [];
else
    % Find scans with tseries files that match those specified in
    % params.tseriesfiles. Use only those scans and the corresponding
    % params.   
    if strcmp(params.tseriesfiles,'any')
        params.tseriesfiles = tseriesfiles;
    end
    match = [];
    scanList = [];
    for s = params.scanList
        m = find(strcmp(tseriesfiles{s},{params.tseriesfiles{:}}));
        if m
            scanList = [scanList,s];
            match = [match,m];
        end
    end
    newparams.scanList = scanList;
    newparams.shiftList = params.shiftList(match);
    newparams.reverseList = params.reverseList(match);
    newparams.tseriesfiles = {params.tseriesfiles{match}};
    newparams.baseScan = params.baseScan;
    newparams.groupName = groupName;
    newparams.aveGroupName = params.aveGroupName;
	newparams.interpMethod = params.interpMethod;
    if ~isfield(params,'description') || isempty(params.description) || isequal(params.description,['Average from ',groupName,' of scans: ',num2str(1:nScans)])
      newparams.description = ['Average from ',groupName,' of scans: ',num2str(scanList)];
    else
      newparams.description = params.description;
    end
    newparams.fileName = [];
end
