% fslFugueTSeriesGUI.m
%
%        $Id: fslFugueTSeriesGUI.m 1833 2010-11-13 18:37:14Z julien $
%      usage: fslFugueTSeriesGUI(varargin)
%         by: julien besle (from motionCompTSeriesGUI
%       date: 11/05/12
%    purpose: gets the parameters for fslFugueTSeries
%
% params = fslFugueTSeriesGUI('groupName','Raw');
% params = fslFugueTSeriesGUI('params',params);
% params = fslFugueTSeriesGUI('groupName','Raw','params',params);
function params = fslFugueTSeriesGUI(varargin)

params = [];
% check arguments
if ~any(nargin == [2 4 6 8])
  help fslFugueTSeriesGUI
  return
end

defaultParams = 0;

% Parse varargin
v = [];
for index = 1:2:length(varargin)
  field = varargin{index};
  val = varargin{index+1};
  switch field
    case 'groupName'
      groupName = val;
    case 'params'
      params = val;
    case 'defaultParams'
      defaultParams = val;
    case 'scanList'
      scanList = val;
    case 'v'
      v = val;
    otherwise
      mrWarnDlg('Invalid initialization argument')
  end
end

% Error if neither params nor groupName specified
if ieNotDefined('params') & ieNotDefined('groupName')
  mrErrorDlg('Must initialize with either params or groupName.');
end

% If groupName not passed then set it according to params.groupName (which
% we now know must exist when groupName does not).
if ieNotDefined('groupName')
  if ~isfield(params,'groupName')
    mrErrorDlg('Must initialize with params.groupName.');
  end
  groupName = params.groupName;
end

% Group names, group number, and nScans
groupNames = viewGet([],'groupNames');
undistortedGroupNames = groupNames;
if ~any(strcmp('MotionComp',undistortedGroupNames))
  undistortedGroupNames{length(undistortedGroupNames)+1} = 'MotionComp';
end
groupNum = viewGet([],'groupNum',groupName);
if isempty(groupNum)
  mrErrorDlg('group ',groupName,' not found.');
end
nScans = viewGet([],'nscans',groupNum);

% Reconcile/initialize params with current status of group and ensure that
% it has the required fields.
if ieNotDefined('params')
  params = fslFugueTSeriesReconcileParams(groupName);
else
  params = fslFugueTSeriesReconcileParams(groupName,params);
end

% Initialize tseriesfiles
tseriesfiles = cell(1,nScans);
for scan = 1:nScans
  tseriesfiles{scan} = viewGet([],'tseriesFile',scan,groupNum);
end

% Initialize descriptions
descriptions = cell(1,nScans);
for scan = 1:nScans
  m = find(scan == params.targetScans);
  if m
    descriptions{scan} = params.descriptions{m};
  else
    descriptions{scan} = ['Motion compensation of ',groupName,' scan ',num2str(scan)];
  end
end

% Initialize group popup
undistortedGroupNames = {'New',undistortedGroupNames{:}};
undistortedGroupName = params.undistortedGroupName;
undistortedGroupNum = find(strcmp(undistortedGroupName,undistortedGroupNames), 1);
if isempty(undistortedGroupNum)
  undistortedGroupNames{length(undistortedGroupNames)+1} = undistortedGroupName;
  undistortedGroupNum = length(undistortedGroupNames);
end

% Initialize base frame popup
unwarpDirectionStrings = {'x','x-','y','y-','z','z-'};

keepAsking=true;
while keepAsking

  % Initialize include (targetScans)
  includeScans = cell(1,nScans);
  for i = 1:length(params.targetScans)
    includeScans{params.targetScans(i)} = 1;
  end

  % set up params
  paramsInfo = {};
  paramsInfo{end+1} = {'undistortedGroupName',putOnTopOfList(undistortedGroupName,undistortedGroupNames),'type=popupmenu','Group name to put the undistorted timeseries scans into'};
  paramsInfo{end+1} = {'fieldMapFile',params.fieldMapFile,'type=string','File containing the B0 field map (in rad/s)'};
  paramsInfo{end+1} = {'chooseFieldMapFile',params.chooseFieldMapFile,'callback',@chooseFieldMapFile,'buttonString=Choose Field Map File','passParams=1','type=pushbutton', 'Click to browse and choose the field map file'};
  paramsInfo{end+1} = {'multiplyBy2pi',params.multiplyBy2pi,'type=checkbox', 'Check this if the field map file is in s-1 (Philips scanner), to transform into rad/s'};
  paramsInfo{end+1} = {'dwellTime',params.dwellTime,'minmax=[0 1]', 'set the EPI dwell time per phase-encode line - same as echo spacing - (in seconds)'};
  paramsInfo{end+1} = {'unwarpDirection',putOnTopOfList(params.unwarpDirection,unwarpDirectionStrings),'type=popupmenu', 'specifies direction of warping (fat shift direction of foldover/phase-encoding axis) e.g. if the phase-encoding axis is RL and the fat shift direction is R, the value should be x-'};
  paramsInfo{end+1} = {'include',includeScans,'type=checkbox','group=scanNum', 'Check this to include a particular scan, uncheck to skip'};
  paramsInfo{end+1} = {'scanNum',1,'incdec=[-1 1]',sprintf('minmax=[1 %i]',nScans),'Scan selector. Use this to choose which scans to include'};
  paramsInfo{end+1} = {'tseriesfiles',tseriesfiles,'group=scanNum','type=String','editable=0','Filename of scan'};
  paramsInfo{end+1} = {'descriptions',descriptions,'group=scanNum','type=String','editable=0','Scan description'};

  % provide parameters for masking if passed in a view that has ROIs loaded
  if isempty(v)
    roiNames = [];
  else
    roiNames = viewGet(v,'roiNames');
  end
  if isempty(roiNames)
    roiMaskOption='enable=0';
    roiNames={'No ROI loaded'};
  else
    roiMaskOption='enable=1';
  end
  paramsInfo{end+1} = {'useMask',0,'type=checkbox',roiMaskOption,'Use a mask for the field map values: specify an ROI and field map outisde this ROI will be set to 0. This is useful to avoid artefacts from large phase gradients just outside the brain.'};
  paramsInfo{end+1} = {'maskROI',roiNames,'The name of the ROI that you want to use as your mask','contingent=useMask'};

  % put up dialog
  if defaultParams
    params = mrParamsDefault(paramsInfo);
    % if a scan list has been specified then
    % only select the scans from the scan list
    if ~ieNotDefined('scanList') 
      if (min(scanList) < 1) || (max(scanList) > nScans)
        mrErrorDlg('(fslFugueTSeriesGUI) scanList include scans that do not exist');
      end
      params.include(:) = 0;
      params.include(scanList) = 1;
    end
  else
    params = mrParamsDialog(paramsInfo,'Set FSL FUGUE parameters');
  end

  if ~isempty(params)
    params.targetScans = find(params.include);
    
    if strcmp(params.fieldMapFile,'empty')
      mrWarnDlg('(fslFugueTSeriesGUI) you must choose a field map file')
    else
      hdr=mlrImageReadNiftiHeader(params.fieldMapFile);
      if any(viewGet(v,'scanDims',params.targetScans(1))~=hdr.dim(2:4)')
        mrWarnDlg('(fslFugueTSeriesGUI) The dimensions of the field map file are not compatible with these scans. Please choose another file.')
      else
        % do some clean up on the parameters
        tseriesfiles = {}; descriptions = {};
        for i = 1:length(params.targetScans)
          tseriesfiles{i} = params.tseriesfiles{params.targetScans(i)};
          descriptions{i} = params.descriptions{params.targetScans(i)};
        end
        params.tseriesfiles = tseriesfiles;
        params.descriptions = descriptions;
        params.groupName=groupName;
        % remove extraneous fields
        params = rmfield(params,'include');
        params = rmfield(params,'paramInfo');
        params = rmfield(params,'scanNum');
        % order fields
        params = orderfields(params);
        keepAsking=false;
      end
    end
  else
    keepAsking=false;
    params=[];
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% callback for setting crop region, gets the correct volume
% and passes it to selectCropRegion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function chooseFieldMapFile = chooseFieldMapFile(params)

if ~isempty(params.fieldMapFile)
  fieldMapPathname = fileparts(params.fieldMapFile);
else
  fieldMapPathname = pwd;
end

[fieldMapFilename fieldMapPathname] = uigetfile({'*.img','*.nii'},'Field map file',fieldMapPathname);
if isnumeric(fieldMapFilename)
  chooseFieldMapFile = false;
  return
else
  newParams.fieldMapFile = [fieldMapPathname fieldMapFilename];
  mrParamsSet(newParams);
  chooseFieldMapFile = 1;
end

