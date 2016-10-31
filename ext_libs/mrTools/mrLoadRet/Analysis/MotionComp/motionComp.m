function [view  params] = motionComp(view,params,varargin)
%
% view = motionComp(view,[params])
%
% Robust 3D rigid body motion compensation. First, performs within scan
% motion compensation on the base scan. Then, computes the mean (over time)
% of the motion-compensated base scan. Finally, motion corrects each
% individual frame of the target scans to the mean motion-corrected base
% scan.
%
% params: Optional initial parameters (default: user is prompted via
%    motionCompGUI). Params must be a structure with all of the
%    following fields.
% targetScans: which scans to apply motion compensation.
%    Default all scans.
% baseFrame: string ('first', last', or 'mean') specifying frame in
%    baseScan to which the rest of the baseScan is registered (ignoring
%    junk frames).
%    Default 'first'.
% robustFlag: use robust M-estimator for motion estimates.
%    Default 0.
% tSmooth: size of the temporal window (in total number of frames)
%    over which to smooth prior to motion estimation
%    Default 0
% driftCorrection: divide each frame by the mean intensity over the volume
%    Default: 1.
% gradIntensityCorrection: intensity gradient correction before registration
%    Default 0.
% crop specifies border size to crop/ignore around all sides of the volume.
%    Should be of the form [ymin xmin zmin; ymax xmax zmax]
%    Default [].
% niters: number of iterations in the motion estimation.
%    Default 10.
% sliceTimeCorrection: Performs slice time correction if True.
% sliceTimeString: Specifies slice time correction:
%    'beginning of TR'
%    'middle of TR'
%    'end of TR' (default)
% interpMethod: 'nearest', 'linear', 'cubic', or 'spline' as in interp3.
% tseriesfiles: cell array of strings the same length as targetScans,
%    specifying the tseries filenames. Or 'any' to allow any match with the
%    tseries files.
%
% If you change this function make parallel changes in:
%   MotionCompBetweenScans, motionCompWithScan
%
%
% Also saves an auxillary file (tseriesfileName.mat) that can be used to
% recompute as follows:
% >> load FileName.mat
% >> eval(evalstr)
%
%
% Examples:
%
% params = motionCompGUImrParams('groupName','Raw');
% view = newView;
% view = motionComp(view,params);
%
% view = motionComp(view);
%
% You can also just get a parameters structure in the following
% ways:
%
%  v = newView;
%  [v params] = motionComp(v,[],'justGetParams=1');
%  [v params] = motionComp(v,[],'justGetParams=1','defaultParams=1');
%  [v params] = motionComp(v,[],'justGetParams=1','defaultParams=1','scanList=[1 2]');
%
%
% DJH 7/2006, updated to MLR4
% jlg 11/2006, save originalFilename and GroupName in scanParams
% Get analysis parameters from motionCompGUI.
% clear time series on each iteration to avoid running out of
% memory
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
  mrWarnDlg('(motionComp) No scans in group');
  return
end

% Get analysis parameters from motionCompGUI.
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
baseScan = params.baseScan;
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

% Error if the number of slices is too small (the derivative computation
% discards the borders in z, 2 slices at the begining and 2 more at the
% end).
baseDataSize = viewGet(viewBase,'datasize',baseScan);
if baseDataSize(3) < 8
  mrErrorDlg('Motion compensation requires at least 8 slices');
end

% Ignore scans if data size is different from that for base scan.
% *** get rid of this by using warpAffine3 below to resample each scan to base scan size
for scanNum = targetScans
  datasize = viewGet(viewBase,'datasize',scanNum);
  if (datasize ~= baseDataSize)
    mrWarnDlg(['(motionComp) Ignoring scan ',num2str(scanNum),'. Motion compensation requires at the datasize to match the base scan']);
    targetScans = targetScans(find(targetScans ~= scanNum));
  end
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Within scan motion compensation on base scan %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load tseries 
scanNum = baseScan;
tseries = loadTSeries(viewBase,scanNum,'all');
junkFrames = viewGet(viewBase,'junkframes',scanNum);
nFrames = viewGet(viewBase,'nFrames',scanNum);
totalFrames = viewGet(viewBase,'totalFrames',scanNum);
if sliceTimeCorrection
  sliceTimes = viewGet(viewBase,'sliceTimes',scanNum);
else
  sliceTimes = [];
end
[mask view] = motionCompGetMask(view,params,scanNum,groupNum);

% Initialize the warped time series to zeros.
% same size as orig tseries
warpedTseries = zeros(size(tseries));

% Preprocess (drift correction, intensit gradient correction, temporal smoothing)
% also correct crop, get slice times, and extract base volume.
[tseriesTemp,crop,newSliceTimes,baseVol,baseF] = motionCompPreprocessing(tseries,params,junkFrames,nFrames,totalFrames,sliceTimes,mask);

% Loop: computing motion estimates and warping the volumes to
% compensate for the motion in each temporal frame.
waitHandle = mrWaitBar(0,['(motionComp) Computing within scan motion compensation for base scan ',num2str(scanNum),'.  Please wait...']);
% Note that the correct transform is just identify here, since we are
% aligning all frames to a base frame from the same run
M = eye(4);
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
        M = estMotionInterp3(baseVol,tseriesTemp,1,frame,niters,Minitial,newSliceTimes,1,robust,0,crop);
      else
        M = estMotionInterp3(tseriesTemp,tseriesTemp,baseF,frame,niters,Minitial,newSliceTimes,1,robust,0,crop);
      end
    else
      M = estMotionIter3(baseVol,tseriesTemp(:,:,:,frame),niters,Minitial,1,robust,0,crop);
    end
  end
  % Warp the volume
  if sliceTimeCorrection
    warpedTseries(:,:,:,frame) = warpAffineInterp3(tseriesTemp,frame,M,newSliceTimes,NaN,interpMethod);
  else
    warpedTseries(:,:,:,frame) = warpAffine3(tseriesTemp(:,:,:,frame),M,NaN,0,interpMethod);
  end
end
mrCloseDlg(waitHandle);

% Finally, compute mean over time (ignoring junkFrames)
baseMean = nanmean(warpedTseries(:,:,:,junkFrames+1:junkFrames+nFrames),4);

% clear temporary tseries
clear tseries tseriesTemp warpedTseries;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through target scans, perform motion estimation for each frame, %
% warp according to motion estimates and save new tseries.             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the base qform and sform.  the warpedTSeries will inheret these
baseQform = viewGet(viewBase, 'scanqform', baseScan, groupName);
baseSform = viewGet(viewBase, 'scansform', baseScan, groupName);

for s = 1:length(targetScans)
  scanNum = targetScans(s);
  
  % Load tseries
  tseries = loadTSeries(viewBase,scanNum,'all');
  junkFrames = viewGet(viewBase,'junkframes',scanNum);
  nFrames = viewGet(viewBase,'nFrames',scanNum);
  totalFrames = viewGet(viewBase,'totalFrames',scanNum); 
  scan2scan = viewGet(viewBase, 'scan2scan', baseScan, groupNum, scanNum, groupNum);
  [mask view] = motionCompGetMask(view,params,scanNum,groupNum);
  if sliceTimeCorrection
    sliceTimes = viewGet(viewBase,'sliceTimes',scanNum);
  else
    sliceTimes = [];
  end
  
  % Initialize the warped time series to zeros. Need to re-initialize this
  % for each scan because number of frames can differ. 
  warpedTseries = zeros(size(tseries));

  % Preprocess (drift correction, intensit gradient correction, temporal smoothing)
  % also correct crop, get slice times, and extract base volume.
  [tseriesTemp,crop,newSliceTimes] = motionCompPreprocessing(tseries,params,junkFrames,nFrames,totalFrames,sliceTimes,mask);
  
  % Loop through frames of target scan and estimate motion params
  waitHandle = mrWaitBar(0,['(motionComp) Computing motion estimates for scan ',num2str(scanNum),'.  Please wait...']);
  % Note that the correct transform is not identity, since we are
  % aligning all frames to a base frame, which could come from a
  % different run
  M = scan2scan;
  transforms = cell(1,totalFrames);
  for frame = 1:totalFrames
    mrWaitBar(frame/totalFrames,waitHandle)
    if (frame <= junkFrames)
      Minitial = scan2scan;
    else
      Minitial = M;
    end
    % Compute rigid-body motion estimate with respect to baseMean
    if sliceTimeCorrection
      M = estMotionInterp3(baseMean,tseriesTemp,1,frame,niters,Minitial,newSliceTimes,1,robust,0,crop);
    else
      M = estMotionIter3(baseMean,tseriesTemp(:,:,:,frame),niters,Minitial,1,robust,0,crop);
    end
    % Collect the transform
    transforms{frame} = M;
  end
  clear tseriesTemp;
  mrCloseDlg(waitHandle);
  
  % warp the images according to the motion estimates
  waitHandle = mrWaitBar(0,['(motionComp) Warping image volumes for scan ',num2str(scanNum),'.  Please wait...']);
  for frame = 1:totalFrames
    mrWaitBar(frame/totalFrames,waitHandle)
    if sliceTimeCorrection
      warpedTseries(:,:,:,frame) = warpAffineInterp3(tseries,frame,transforms{frame},newSliceTimes,NaN,interpMethod);
    else
      warpedTseries(:,:,:,frame) = warpAffine3(tseries(:,:,:,frame),transforms{frame},NaN,0,interpMethod);
    end
  end
  mrCloseDlg(waitHandle);
  clear tseries;
  
  % Save tseries with modified nifti header
  scanParams = viewGet(viewBase,'scanParams',scanNum);
  scanParams.description = ['Full ' descriptions{s}];
  scanParams.fileName = [];
  scanParams.originalFileName{1} = viewGet(viewBase,'tseriesfile',scanNum);
  scanParams.originalGroupName{1} = viewGet(viewBase,'groupName');

  % set the qform and sform to that of the baseScan
  scanParams.niftiHdr = cbiSetNiftiQform(scanParams.niftiHdr, baseQform);
  scanParams.niftiHdr = cbiSetNiftiSform(scanParams.niftiHdr, baseSform);
  
  [viewMotionComp,tseriesFileName] = saveNewTSeries(viewMotionComp,warpedTseries,scanParams,scanParams.niftiHdr);

  % Save evalstring for recomputing and params
  evalstr = ['view = newView(','''','Volume','''','); view = motionComp(view,params);'];
  [pathstr,filename] = fileparts(tseriesFileName);
  tseriesdir = viewGet(viewMotionComp,'tseriesdir');
  save(fullfile(tseriesdir,filename),'evalstr','params','transforms','tseriesFileName');

  % clear temporary tseries
  clear warpedTseries;

end

% Delete temporary views
deleteView(viewBase);
deleteView(viewMotionComp);

return

