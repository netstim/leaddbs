function view = sliceTimeCorrect(view,params)
%
% view = sliceTimeCorrect(view,[params])
%
% Slice time correction without motion compensation
%
% params: Optional initial parameters (default: user is prompted via
%    motionCompGUI). Params must be a structure with all of the
%    following fields.
% targetScans: which scans to apply motion compensation.
%    Default all scans.
% sliceTimeCorrection: must be True
% interpMethod: 'nearest', 'linear', 'cubic', or 'spline' as in interp3.
% tseriesfiles: cell array of strings the same length as targetScans,
%    specifying the tseries filenames. Or 'any' to allow any match with the
%    tseries files.
%
% If you change this function make parallel changes in:
%   motionCompBetweenScans, motionComp, motionCompWithinScan
%
%
% Examples:
%
% params = motionCompGUI('groupName','Raw');
% view = newView;
% view = sliceTimeCorrect(view,params);
%
% view = sliceTimeCorrect(view);
%
%
% DJH 7/2007, modified from motionCompWithinScan

% Get analysis parameters from motionCompGUI.
nScans = viewGet(view,'nScans');
if (nScans == 0)
  mrWarnDlg('(sliceTimeCorrect) No scans in group');
  return
end

if ieNotDefined('params')
  % Initialize analysis parameters with default values
  params = motionCompGUI('groupName',viewGet(view,'groupName'));
else
  % Reconcile params with current status of group and ensure that it has
  % the required fields.
  params = motionCompReconcileParams(params.groupName,params);
end

% Abort if params empty or if sliceTimeCorrection not selected
if ieNotDefined('params')
  mrMsgBox('Slice time correction cancelled');
  return
end
if ~params.sliceTimeCorrection
  mrMsgBox('Slice time correction cancelled');
  return
end

% Retrieve parameters
targetScans = params.targetScans;
interpMethod = params.interpMethod;
sliceTimeString = params.sliceTimeString;
groupName = params.groupName;
motionCompGroupName = params.motionCompGroupName;
descriptions  = params.descriptions;
tseriesfiles = params.tseriesfiles;

% Open new view with the base group
viewBase = newView(viewGet(view,'viewType'));
groupNum = viewGet(viewBase,'groupNum',groupName);
if (groupNum == 0)
  mrErrorDlg('motionComp: ',groupName,' does not exist.');
end
viewBase = viewSet(viewBase,'currentGroup',groupNum);

% Open new view and set its group to the motion comp group name. Create the
% group if necessary.
viewMotionComp = newView(viewGet(view,'viewType'));
motionCompGroupNum = viewGet(viewMotionComp,'groupNum',motionCompGroupName);
if isempty(motionCompGroupNum)
  view = viewSet(view,'newgroup',motionCompGroupName);
  motionCompGroupNum = viewGet(viewMotionComp,'groupNum',motionCompGroupName);
end
viewMotionComp = viewSet(viewMotionComp,'currentGroup',motionCompGroupNum);

% Loop through target scans, performing slice time correction for each
% frame, and save new tseries.
for s = 1:length(targetScans)
  scanNum = targetScans(s);

  % Load tseries 
  tseries = loadTSeries(viewBase,scanNum,'all');
  junkFrames = viewGet(viewBase,'junkframes',scanNum);
  nFrames = viewGet(viewBase,'nFrames',scanNum);
  totalFrames = viewGet(viewBase,'totalFrames',scanNum);
  
  % Get slice times and replicate the last frame for slice time correction
  warpedTseries = zeros(size(tseries));
  
  % Replicate last frame of tseries
  tseries(:,:,:,end+1) = tseries(:,:,:,end);
  
  % Get slice times
  sliceTimes = viewGet(viewBase,'sliceTimes',scanNum);
  switch sliceTimeString
    case 'end of TR'
      sliceTimes = sliceTimes;
    case 'middle of TR'
      sliceTimes = sliceTimes - 0.5;
    case 'beginning of TR'
      sliceTimes = sliceTimes - 1;
    otherwise
      mrErrorDlg('Invalid slice times');
  end

  % Make grids for interpn
  [Ny Nx Nz Nt] = size(tseries);
  [ygrid,xgrid,zgrid,tgrid] = ndgrid(1:Ny,1:Nx,1:Nz,0);
  % Slice time correction interpolated to the end of each TR. Would have to
  % add 1/2 to interpolate to the middle of each TR but that busts the last
  % frame of the tseries.
  for slice = 1:size(tgrid,3)
    tgrid(:,:,slice) = tgrid(:,:,slice) - sliceTimes(slice);
  end

  % Loop through frames
  waitHandle = mrWaitBar(0,['Computing slice time correction for scan ',num2str(scanNum),'.  Please wait...']);
  for frame = 1:totalFrames
    mrWaitBar(frame/totalFrames,waitHandle)
    tgrid = tgrid + 1;
    warpedTseries(:,:,:,frame) = interpn(tseries,ygrid,xgrid,zgrid,tgrid,interpMethod,NaN);
  end
  mrCloseDlg(waitHandle);

  % Save tseries with modified nifti header
  scanParams = viewGet(viewBase,'scanParams',scanNum);
  scanParams.description = ['Within ' descriptions{s}];
  scanParams.fileName = [];
  scanParams.originalFileName{1} = viewGet(viewBase,'tseriesfile',scanNum);
  scanParams.originalGroupName{1} = viewGet(viewBase,'groupName');
  [viewMotionComp,tseriesFileName] = saveNewTSeries(viewMotionComp,warpedTseries,scanParams,scanParams.niftiHdr);

  % Save evalstring for recomputing and params
  evalstr = ['view = newView(','''','Volume','''','); view = sliceTimeCorrect(view,params);'];
  [pathstr,filename,ext] = fileparts(tseriesFileName);
  tseriesdir = viewGet(viewMotionComp,'tseriesdir');
  save(fullfile(tseriesdir,filename),'evalstr','params','tseriesFileName');
  
  % clear temporary tseries
  clear tseries tseriesIC warpedTseries tseriesCrop;

end

% Delete temporary views
deleteView(viewBase);
deleteView(viewMotionComp);

return
