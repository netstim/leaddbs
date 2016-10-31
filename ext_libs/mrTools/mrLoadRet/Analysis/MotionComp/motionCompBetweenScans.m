function [view params] = motionCompBetweenScans(view,params,varargin)
%
% view = motionCompBetweenScans(view,[params],[justGetParams])
%
% Robust 3D rigid body motion compensation between different scans.
% Computes the mean over time of each scan, and then computes the
% registration between those mean volumes. Finally, rather than
% interpolating the time series, this function simply copies the tseries
% file, overwriting the sform of the nifti header.
%
% params: Optional initial parameters (default: user is prompted via
%    motionCompGUI). Params must be a structure.  See motionComp.m for
%    details.
% 
% justGetParams: option to return param structure without running the
%    motion compensation (see motionComp.m for details and examples).
%
% If you change this function make parallel changes in:
%    motionComp, motionCompWithinScan
%
% Also saves an auxillary file (tseriesfileName.mat) that can be used to
% recompute as follows:
% >> load FileName.mat
% >> eval(evalstr)
%
%
% Examples:
%
% params = motionCompGUI('groupName','Raw');
% view = newView;
% view = motionCompBetweenScans(view,params);
%
% view = motionCompBetweenScans(view);
%
%
% DJH 7/2006, updated to MLR4
% jlg 11/2006, save originalFilename and GroupName in scanParams
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
    mrWarnDlg('(motionCompBetweenScans) No scans in group');
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
for scanNum = targetScans
  datasize = viewGet(viewBase,'datasize',scanNum);
  if (datasize ~= baseDataSize)
    mrWarnDlg(['Ignoring scan ',num2str(scanNum),'. Motion compensation requires at least 8 slices']);
    targetScans = targetScans(find(targetScans ~= scanNum));
  end
end

% Open new view and set its group to the motion comp group name. Create the
% group if necessary.
viewMotionComp = newView(viewGet(view,'viewType'));
motionCompGroupNum = viewGet(viewMotionComp,'groupNum',motionCompGroupName);
if isempty(motionCompGroupNum)
  view = viewSet(view,'newgroup',motionCompGroupName);
  motionCompGroupNum = viewGet(viewMotionComp,'groupNum',motionCompGroupName);
end
viewMotionComp = viewSet(viewMotionComp,'currentGroup',motionCompGroupNum);

% Get base scanParams
scanParams = viewGet(viewBase,'scanParams',baseScan);
%baseSform = scanParams.niftiHdr.sform44;

% Open wait bar
waitHandle = mrWaitBar(0,'Computing between scan motion compensation.  Please wait...');

% Compute baseMean
% Load tseries and dump junk frames
tseries = loadTSeries(viewBase,baseScan,'all');
junkFrames = viewGet(viewBase,'junkframes',baseScan);
nFrames = viewGet(viewBase,'nFrames',baseScan);
% Compute mean of base scan
baseMean = nanmean(tseries(:,:,:,junkFrames+1:junkFrames+nFrames),4);
[mask view] = motionCompGetMask(view,params,baseScan,groupNum);

% Preprocess (drift correction, intensit gradient correction, temporal smoothing)
% also correct crop, get slice times, and extract base volume.
[baseMean,crop] = motionCompPreprocessing(baseMean,params,0,1,1,[],mask);

% Loop through target scans, compute motion estimate transform, warp
% according to the transform, and save new tseries.
for s = 1:length(targetScans)
  scanNum = targetScans(s);
  
  % Update wait bar
  mrWaitBar(s/length(targetScans),waitHandle);
  
  % Load tseries 
  tseries = loadTSeries(viewBase,scanNum,'all');
  
  % Load the transform that brings the scan into rough alignment with
  % the base
  scan2scan = viewGet(viewBase, 'scan2scan', baseScan, groupNum, scanNum, groupNum);
  swapXY = [0 1 0 0;1 0 0 0;0 0 1 0; 0 0 0 1];
  Minitial = swapXY*scan2scan*swapXY;
  
  if (scanNum == baseScan)
    % Do not bother computing registration for baseScan (just copy it)
    transform = eye(4);
    warpedTseries = tseries;
  else
    junkFrames = viewGet(viewBase,'junkframes',scanNum);
    nFrames = viewGet(viewBase,'nFrames',scanNum);
    totalFrames = viewGet(viewBase,'totalFrames',scanNum);

    % Initialize the warped time series to zeros. Need to re-initialize this
    % for each scan because number of frames can differ. 
    warpedTseries = zeros(size(tseries));
    
    % get mask
    [mask view] = motionCompGetMask(view,params,scanNum,groupNum);

    % Compute mean of target scan after dumping junk frames
    targetMean = nanmean(tseries(:,:,:,junkFrames+1:junkFrames+nFrames),4);
    % Preprocess (drift correction, intensit gradient correction, temporal smoothing)
    % also correct crop, get slice times, and extract base volume.
    targetMean = motionCompPreprocessing(targetMean,params,0,1,1,[],mask);

    % Estimate motion between mean volumes
    disp(sprintf('Computing motion estimates for scan %i ...',scanNum))
    %transform = estMotionIter3(baseMean,targetMean,niters,eye(4),1,robust,0,crop);
    transform = estMotionIter3(baseMean,targetMean,niters,Minitial,1,robust,0,crop);
    disp(sprintf('Computing motion estimates for scan %i ...done',scanNum))
    
    % jg: warp/transform each frame of the time series, rather than saving
    % the transformation in the sform
    disp(['Warping images for scan ',num2str(scanNum)])
    for frame = 1:totalFrames
      warpedTseries(:,:,:,frame) = warpAffine3(tseries(:,:,:,frame),transform,NaN,0,interpMethod);
    end
    disp(['Warping images for scan ',num2str(scanNum),'... done'])
  end
    
  % Save tseries with modified nifti header
  scanParams = viewGet(viewBase,'scanParams',scanNum);
  scanParams.description = ['Between ' descriptions{s}];
  scanParams.junkFrames = 0;
  scanParams.nFrames = nFrames;
  scanParams.fileName = [];

  % get the base scan's qform and sform -- the warpedTSeries will inheret these
  baseQform = viewGet(viewBase, 'scanqform', baseScan, groupName);
  baseSform = viewGet(viewBase, 'scansform', baseScan, groupName);

  % set the qform and sform to that of the baseScan
  scanParams.niftiHdr = cbiSetNiftiQform(scanParams.niftiHdr, baseQform);
  scanParams.niftiHdr = cbiSetNiftiSform(scanParams.niftiHdr, baseSform);

  
  scanParams.originalFileName{1} = viewGet(viewBase,'tseriesfile',scanNum);
  scanParams.originalGroupName{1} = viewGet(viewBase,'groupName');
  [viewMotionComp,tseriesFileName] = saveNewTSeries(viewMotionComp,warpedTseries,scanParams,scanParams.niftiHdr);
  
  % Save evalstring for recomputing and params
  evalstr = ['view = newView(','''','Volume','''','); view = motionCompBetweenScans(view,params);'];
  [pathstr,filename,ext] = fileparts(tseriesFileName);
  tseriesdir = viewGet(viewMotionComp,'tseriesdir');
  save(fullfile(tseriesdir,filename),'evalstr','params','transform','tseriesFileName');
  
  % clear temporary tseries
  clear tseries  warpedTseries;

end 
% Close waitbar
mrCloseDlg(waitHandle);

% Delete temporary views
deleteView(viewBase);
deleteView(viewMotionComp);

return
