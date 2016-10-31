function newparams = fslFugueTSeriesReconcileParams(groupName,params,varargin)
% params = fslFugueTSeriesReconcileParams(groupName,[params])
%
% Checks for consistency between fslFugueTSeries parameters and current
% status of group (nscans and tseries filenames). Reconciles params by
% matching tseries filenames.
%
% params: Optional initial parameters (see fslFugueTSeries).
% default if not passed: initialize new params with default values.
%
% jb 11/05/2012 (from motionCompReconcileParams)

groupNum = viewGet([],'groupNum',groupName);
nScans = viewGet([],'nscans',groupNum);

% Get tseries filenames for this group
tseriesfiles = cell(1,nScans);
for scan = 1:nScans
  tseriesfiles{scan} = viewGet([],'tseriesFile',scan,groupNum);
  oldDescription = viewGet([],'description',scan,groupNum);
  descriptions{scan} = ['B0 correction of ',groupName,' scan ',num2str(scan) ': ' oldDescription];
end

if ieNotDefined('params')
  % Use default params
  newparams.fieldMapFile = 'empty';
  newparams.chooseFieldMapFile = 0;
  newparams.multiplyBy2pi = 1;
  newparams.groupName = groupName;
  newparams.dwellTime = 0.000477;
  newparams.unwarpDirection = 'x';
  newparams.undistortedGroupName = 'B0 corrected';
  newparams.targetScans = 1:nScans;
  newparams.tseriesfiles = tseriesfiles;
  newparams.descriptions = descriptions;
else
  % Set newparams according to params, reconciling with tseries files.
  newparams.fieldMapFile = params.fieldMapFile;
  newparams.chooseFieldMapFile = params.chooseFieldMapFile;
  newparams.multiplyBy2pi = params.multiplyBy2pi;
  newparams.groupName = params.groupName;
  newparams.dwellTime = params.dwellTime;
  newparams.unwarpDirection = params.unwarpDirection;
  newparams.undistortedGroupName = params.undistortedGroupName;
  % Find scans with tseries files that match those specified in
  % params.tseriesfiles. Use only those scans and the corresponding
  % params.
  if strcmp(params.tseriesfiles,'any')
    params.tseriesfiles = tseriesfiles;
    params.descriptions = descriptions;
  end
  match = [];
  targetScans = [];
  for s = params.targetScans
    m = find(strcmp(tseriesfiles{s},{params.tseriesfiles{:}}));
    if m
      targetScans = [targetScans,s];
      match = [match,m];
    end
  end
  newparams.targetScans = targetScans;
  newparams.tseriesfiles = {params.tseriesfiles{match}};
  newparams.descriptions = {params.descriptions{match}};
end
