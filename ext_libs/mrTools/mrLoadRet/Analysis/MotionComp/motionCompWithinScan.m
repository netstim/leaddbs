function [view params] = motionCompWithinScan(view,params,varargin)
%
% view = motionCompWthinScan(view,[params],[justGetParams])
%
% Robust 3D rigid body motion compensation across frames within each
%    scan.
%
% params: Optional initial parameters (default: user is prompted via
%    motionCompGUI). Params must be a structure.  See motionComp.m for
%    details.
% 
% justGetParams: option to return param structure without running the
%    motion compensation (see motionComp.m for details and examples).
%
% If you change this function make parallel changes in:
%   MotionCompBetweenScans, motionComp
%
% Also saves an auxillary file (tseriesfileName.mat) that can be used to
% recompute as follows:
% >> load FileName.mat
% >> eval(evalstr)
%
% Examples:
%
% params = motionCompGUImrParams('groupName','Raw');
% view = newView;
% view = motionCompWithinScan(view,params);
%
% view = motionCompWithinScan(view);
%
% DJH 7/2006, updated to MLR4
% jlg 11/2006, save originalFilename and GroupName in scanParams,
% clear time series on each iteration to avoid running out of memory
%
% $Id$	
%

% other arguments
eval(evalargs(varargin));
if ieNotDefined('justGetParams'),justGetParams = 0;end
if ieNotDefined('defaultParams'),defaultParams = 0;end
if ieNotDefined('scanList'),scanList = [];end

nScans = viewGet(view,'nScans');
if (nScans == 0)
  mrWarnDlg('(motionCompWithinScans) No scans in group');
  return
end

% Get analysis parameters from motionCompGUImrParams.
if ieNotDefined('params')
  % Initialize analysis parameters with default values
  params = motionCompGUImrParams('groupName',viewGet(view,'groupName'),'defaultParams',defaultParams,'scanList',scanList,'v',view);
else
  % Reconcile params with current status of group and ensure that it has
  % the required fields.
  params = motionCompReconcileParams(params.groupName,params);
end

% Abort if params empty
if ieNotDefined('params')
  mrMsgBox('motion compensation cancelled');
  return
end

% if just getting params then return
if justGetParams,return,end

% Retrieve parameters
baseFrame = params.baseFrame;
targetScans = params.targetScans;
sliceTimeCorrection = params.sliceTimeCorrection;
sliceTimeString = params.sliceTimeString;
robust = params.robust;
niters = params.niters;
interpMethod = params.interpMethod;
groupName = params.groupName;
motionCompGroupName = params.motionCompGroupName;
descriptions  = params.descriptions;
tseriesfiles = params.tseriesfiles;

% Open new view with the base group
viewBase = newView;
groupNum = viewGet(viewBase,'groupNum',groupName);
if (groupNum == 0)
  mrErrorDlg('motionComp: ',groupName,' does not exist.');
end
viewBase = viewSet(viewBase,'currentGroup',groupNum);

% Ignore scans if the number of slices is too small (the derivative
% computation discards the borders in z, 2 slices at the begining and 2
% more at the end).
for scanNum = targetScans
  datasize = viewGet(viewBase,'datasize',scanNum);
  if (datasize(3) < 8)
    mrWarnDlg(['Ignoring scan ',num2str(scanNum),'. Motion compensation requires at least 8 slices']);
    targetScans = targetScans(find(targetScans ~= scanNum));
  end
end
if length(targetScans) == 0
  mrErrorDlg('(motionComWithinScan) No target scans with 8 or more slices');
end

% Open new view and set its group to the motion comp group name. Create the
% group if necessary.
viewMotionComp = newView;
motionCompGroupNum = viewGet(viewMotionComp,'groupNum',motionCompGroupName);
if isempty(motionCompGroupNum)
  view = viewSet(view,'newgroup',motionCompGroupName);
  motionCompGroupNum = viewGet(viewMotionComp,'groupNum',motionCompGroupName);
end
viewMotionComp = viewSet(viewMotionComp,'currentGroup',motionCompGroupNum);

% Loop through target scans, perform motion estimation for each frame, warp
% according to motion estimates and save new tseries.
for s = 1:length(targetScans)
  scanNum = targetScans(s);

  % Load tseries 
  tseries = loadTSeries(viewBase,scanNum,'all');
  junkFrames = viewGet(viewBase,'junkframes',scanNum);
  nFrames = viewGet(viewBase,'nFrames',scanNum);
  totalFrames = viewGet(viewBase,'totalFrames',scanNum);
  [mask view] = motionCompGetMask(view,params,scanNum,groupNum);
  if sliceTimeCorrection
    sliceTimes = viewGet(viewBase,'sliceTimes',scanNum);
  else
    sliceTimes = [];
  end
  
  % Initialize the warped time series to zeros.
  % same size as orig tseries
  warpedTseries = zeros(size(tseries));

  % Preprocess (drift correction, intensit gradient correction, temporal smoothing)
  % also correct crop, get slice times, and extract base volume.
  [tseriesTemp,crop,sliceTimes,baseVol,baseF] = motionCompPreprocessing(tseries,params,junkFrames,nFrames,totalFrames,sliceTimes,mask);

  % Loop: computing motion estimates and warping the volumes to
  % compensate for the motion in each temporal frame.
  waitHandle = mrWaitBar(0,['Computing motion compensation for scan ',num2str(scanNum),'.  Please wait...']);
  M = eye(4);
  transforms = cell(1,totalFrames);
  for frame = 1:totalFrames
    mrWaitBar(frame/totalFrames,waitHandle)
    if (frame == baseF)
      M = eye(4);
    else
      if (frame <= junkFrames)
        Minitial = eye(4);
      else
        Minitial = M;
      end
      % Compute rigid-body motion estimate
      if sliceTimeCorrection
        if strcmp(baseFrame,'mean')
          M = estMotionInterp3(baseVol,tseriesTemp,1,frame,niters,Minitial,sliceTimes,1,robust,0,crop);
        else
          M = estMotionInterp3(tseriesTemp,tseriesTemp,baseF,frame,niters,Minitial,sliceTimes,1,robust,0,crop);
        end
      else
        M = estMotionIter3(baseVol,tseriesTemp(:,:,:,frame),niters,Minitial,1,robust,0,crop);
      end
    end
    % Collect the transform
    transforms{frame} = M;
    % Warp the volume
    if sliceTimeCorrection
      warpedTseries(:,:,:,frame) = warpAffineInterp3(tseries,frame,M,sliceTimes,NaN,interpMethod);
    else
      warpedTseries(:,:,:,frame) = warpAffine3(tseries(:,:,:,frame),M,NaN,0,interpMethod);
    end
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
  evalstr = ['view = newView(','''','Volume','''','); view = motionCompWithinScan(view,params);'];
  [pathstr,filename,ext] = fileparts(tseriesFileName);
  tseriesdir = viewGet(viewMotionComp,'tseriesdir');
  save(fullfile(tseriesdir,filename),'evalstr','params','transforms','tseriesFileName');
  
  % clear temporary tseries
  clear tseries tseriesTemp warpedTseries;

end

% Delete temporary views
deleteView(viewBase);
deleteView(viewMotionComp);

return
