% motionCompGUImrParams.m
%
%        $Id$
%      usage: motionCompGUImrParams(varargin)
%         by: justin gardner
%       date: 08/03/07
%    purpose: replace the old GUI with an mrParams dialog,
%             initial portion that parses input arguments taken
%             from motionComp.m
%
% params = motionCompGUImrParams('groupName','Raw');
% params = motionCompGUImrParams('params',params);
% params = motionCompGUImrParams('groupName','Raw','params',params);
function mrParams = motionCompGUImrParams(varargin)

mrParams = [];
% check arguments
if ~any(nargin == [2 4 6 8])
  help motionCompGUImrParams
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
if ~any(strcmp('MotionComp',groupNames))
  groupNames{length(groupNames)+1} = 'MotionComp';
end
groupNum = viewGet([],'groupNum',groupName);
if isempty(groupNum)
  mrErrorDlg('group ',groupName,' not found.');
end
nScans = viewGet([],'nscans',groupNum);

% Reconcile/initialize params with current status of group and ensure that
% it has the required fields.
if ieNotDefined('params')
  params = motionCompReconcileParams(groupName);
else
  params = motionCompReconcileParams(groupName,params);
end

% Initialize include (targetScans)
includeScans = cell(1,nScans);
for i = 1:length(params.targetScans)
  includeScans{params.targetScans(i)} = 1;
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
makeNewGroupStr = 'Make a new group';
groupNames = {'Make a new group',groupNames{:}};
motionCompGroupName = params.motionCompGroupName;
motionCompGroupNum = find(strcmp(motionCompGroupName,groupNames));
if isempty(motionCompGroupNum)
  groupNames{length(groupNames)+1} = motionCompGroupName;
  motionCompGroupNum = length(groupNames);
end

% Initialize base frame popup
baseFrameStrings = {'first','last','mean'};
% Initialize interp method popup
interpMethodStrings = {'nearest','linear','cubic','spline'};
% Initialize sliceTime strings
sliceTimeStrings = {'beginning of TR','middle of TR','end of TR'};

% set up params
paramsInfo = {};
paramsInfo{1} = {'motionCompGroupName',putOnTopOfList(motionCompGroupName,groupNames),'type=popupmenu','Group name to put the motion compensated scans into'};
paramsInfo{2} = {'interpMethod',putOnTopOfList(params.interpMethod,interpMethodStrings),'type=popupmenu', 'Interpolation method used to warp the images'};
paramsInfo{3} = {'baseScan',params.baseScan,'incdec=[-1 1]',sprintf('minmax=[1 %i]',nScans), 'Specifies the scan that everything else will be aligned to'};
paramsInfo{4} = {'baseFrame',putOnTopOfList(params.baseFrame,baseFrameStrings),'type=popupmenu', 'Specifies which frame (or mean) to align to.  First frame can be bad choice b/c of different T2 contrast'};
paramsInfo{5} = {'niters',params.niters,'incdec=[-1 1]','minmax=[0 inf]', 'How many iterations to estimate the optimal transform (use more for high-res images)'};
paramsInfo{6} = {'sliceTimeCorrection',params.sliceTimeCorrection,'type=checkbox', 'Apply slice time correction along with motion compensation.  Not appropriate for 3D scans.'};
paramsInfo{7} = {'sliceTimeString',putOnTopOfList(params.sliceTimeString,sliceTimeStrings),'type=popupmenu','Which point in time the slices should be aligned to. May loose first and last frames'};
paramsInfo{8} = {'driftCorrection', params.driftCorrection, 'type=checkbox', 'Correction for fluctuations in mean intensity over time. This divides each frame by the average over the (cropped) volume. Note that this is only used to create a temporary time series for better estimation of motion parameters. The actual motion comp volume will be created from the original time series without drift correction applied.'};
paramsInfo{9} = {'tSmooth', params.tSmooth, 'incdec=[-1 1]', 'minmax=[0 10]','round=1' 'Size of boxcar window for temporal smoothing. That is, this averages across tSmooth frames before and after every frame. So, if set to 1 it will create a window of 1 volume on either side of every frame and average 3 frames together. This can be useful to get rid of sudden changes in the image which may not be real motion but are image artifacts (assuming that real head motion is slow).  Like driftCorrection, only applied to estimate head motion, not to final time series'};
paramsInfo{10} = {'gradIntensityCorrection',params.gradIntensityCorrection,'type=checkbox', 'Compensation for intensity gradient before motion estimation. Like driftCorrection this is only applied to estimate head motion, not to final time series. This uses the crop region below'};
paramsInfo{11} = {'crop',params.crop,'callback',@thisSelectCropRegion,'buttonString=Set crop region','passParams=1','type=pushbutton','contingent=gradIntensityCorrection','Crop the images.  This affects the gradIntensityCorrection.  Important for high-res images. The selected crop region is used to sample the intensities and should be draw in a small region of the brain. This does not mean that only the crop region will be used for motion comp (the whole volume will be used). If you want to crop so that only some small region of the volume is used for motion comp you should use the mask option which will be visible if you have an ROI loaded.'};
paramsInfo{12} = {'robust',params.robust,'type=checkbox','contingent=gradIntensityCorrection','Robust contrast estimator, should be used if images are noisy with lots of outliers. This applies to gradIntensityCorrection'};
paramsInfo{13} = {'include',includeScans,'type=checkbox','group=scanNum', 'Check this to include a particular scan, uncheck to skip'};
paramsInfo{14} = {'scanNum',1,'incdec=[-1 1]',sprintf('minmax=[1 %i]',nScans),'Scan selector. Use this to choose which scans to include'};
paramsInfo{15} = {'tseriesfiles',tseriesfiles,'group=scanNum','type=String','editable=0','Filename of scan'};
paramsInfo{16} = {'descriptions',descriptions,'group=scanNum','type=String','editable=0','Scan description'};

% provide parameters for masking if passed in a view that has ROIs loaded
if ~isempty(v)
  roiNames = viewGet(v,'roiNames');
  if ~isempty(roiNames)
    paramsInfo{end+1} = {'useMask',0,'type=checkbox','Use a mask. This is different from the crop region before in that it sets the region which you want to do the motion comp over. That is, you can specify an ROI and then only do motion compensation based on that part of the image that falls within the ROI. THis is useful if you have noise in your image outside the brain.'};
    paramsInfo{end+1} = {'maskROI',roiNames,'The name of the ROI that you want to use as your mask','contingent=useMask'};
% following not implemented yet
%    maskTypes = {'Mask only during motion comp','Remove masked points from tSeries'};
    %    paramsInfo{end+1} = {'maskType',maskTypes,'You can either just compute the motion compensation based on the mask region, but in the final time series keep the points in the image outside the mask, or you can compeltely mask out points of the image in the tSeries in which case they will show up in the final tSeries as nan','contingent=useMask'};
  end
end
% put up dialog
if defaultParams
  mrParams = mrParamsDefault(paramsInfo);
  % if a scan list has been specified then
  % only select the scans from the scan list
  if ~ieNotDefined('scanList') 
    if (min(scanList) < 1) || (max(scanList) > nScans)
      mrErrorDlg('(motionCompGUImrParams) scanList include scans that do not exist');
    end
    mrParams.include(:) = 0;
    mrParams.include(scanList) = 1;
  end
else
  mrParams = mrParamsDialog(paramsInfo,'Set motionComp parameters');
end

if ~isempty(mrParams)
  % save these params to come up as defaults the next time
  mrSetPref('motionCompDefaultParams',mrParams);
  % do some clean up on the parameters
  % first check for a new group
  if isfield(mrParams,'motionCompGroupName') && ~isempty(mrParams.motionCompGroupName) && strcmp(mrParams.motionCompGroupName,makeNewGroupStr)
    paramsInfoNewGroupName = {{'groupName','','Name for new group'}};
    groupName = mrParamsDialog(paramsInfoNewGroupName,'Name for new group');
    if isempty(groupName) mrParams = [];return, end
    if isempty(groupName.groupName) || ~isstr(groupName.groupName)
      disp(sprintf('(motionCompGUImrParams) Invalid groupName'));
      mrParams = [];
      return;
    else
      mrParams.motionCompGroupName = fixBadChars(groupName.groupName);
    end
  end
  % next makes sure it return exactly what the old GUI returned
  mrParams.groupName = params.groupName;
  if isempty(mrParams.crop)
    mrParams.crop = [];
  end
  mrParams.targetScans = find(mrParams.include);
  tseriesfiles = {}; descriptions = {};
  for i = 1:length(mrParams.targetScans)
    tseriesfiles{i} = mrParams.tseriesfiles{mrParams.targetScans(i)};
    descriptions{i} = mrParams.descriptions{mrParams.targetScans(i)};
  end
  mrParams.tseriesfiles = tseriesfiles;
  mrParams.descriptions = descriptions;
  % remove extraneous fields
  mrParams = rmfield(mrParams,'include');
  mrParams = rmfield(mrParams,'paramInfo');
  mrParams = rmfield(mrParams,'scanNum');
  % order fields
  mrParams = orderfields(mrParams);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% callback for setting crop region, gets the correct volume
% and passes it to selectCropRegion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function crop = thisSelectCropRegion(params)

view = newView;
% Load first frame of base scan
volume = loadTSeries(view,params.baseScan,'all',1);
% and get crop region
crop = selectCropRegion(volume);
