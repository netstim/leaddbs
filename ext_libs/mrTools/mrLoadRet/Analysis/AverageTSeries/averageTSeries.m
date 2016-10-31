function [view params] = averageTSeries(view,params,varargin)
%
% function averageTSeries(view,[params])
%
% Computes average of the time series. Time reverses and shifts as
% specified. Output is a new tSeries file added to the specified group, and
% that is in the coordinate frame of the base scan.
%
% Order of operations is: 1) dump junk frames, 2) time shift, 3) time
% reverse, 4) average (transforming to coordinate frame of the base scan).
%
% params: Optional initial parameters (default: user is prompted via
%    averageTSeriesGUI). Params must be a structure with all of the
%    following fields.
% scanList: vector of scan numbers to include in the average.
%    Default: all of the scans.
% tseriesfiles: cell array of strings the same length as scanList,
%    specifying the tseries filenames. Or 'any' to allow any match with the
%    tseries files.
% shiftList: vector of integers the same length as scanList, each entry of
%    which specifies the time shift (in frames, i.e., 1 means a shift of 1
%    TR) for each scan. 
%    Default: all zero (no shift).
% reverseList: vector of the same length as scanList with non-zero entries
%    indicating which scans to time reverse.
%    Default: all zero (no time reversal).
% baseScan: number specifying scan which is used to specify the coordinate
%    frame for the result.
%    Default: first scan
% interpMethod: 'nearest', 'linear', 'cubic', or 'spline' as in interp3.
% groupName: group of scans that will be averaged.
%    Default: current group of view.
% aveGroupName: the new average scan is added to the specified group. If
%    the group specified by aveGroupName does not exist, then a new group
%    is created with this name.
%    Default: 'Averages'.
% description: description string for the new average scan
%    Default: 'Average from <groupName>' of scans: <scanList>'
% fileName: the new tseries is saved in the specified filename.
%    Default: 'tseries-mmddyy-hhmmss.img' or 'tseries-mmddyy-hhmmss.nii'
%    where mmddyy = date and hhmmss = time. Chooses between .img and .nii,
%    based on 'niftiFileExtension' preference.
%
% Also saves an auxillary file (tseriesfileName.mat) that can be used to
% recompute as follows:
% >> load FileName.mat
% >> eval(evalstr)
%
%
% Examples:
%
% params = averageTSeriesGUI('groupName','Raw');
% view = newView;
% view = averageTSeries(view,params);
%
% view = averageTSeries(view);
%
% To just get parameters
% [view params] = averageTSeries(view,[],'justGetParams=1','defaultParams=1','scanList=[1 8]');
%
% djh, 7/2006
% $Id$	

% check arguments
if ~any(nargin == [1 2 3 4 5 6 7 8])
  help averageTSeries
  return
end

% other arguments
eval(evalargs(varargin));
if ieNotDefined('justGetParams'),justGetParams = 0;end
if ieNotDefined('defaultParams'),defaultParams = 0;end

% if we are just getting default parameters then
% get them by calling the reconcile function
if defaultParams
  if ieNotDefined('scanList')
    params = averageTSeriesReconcileParams(viewGet(view,'groupName'),[]);
  else
    params = averageTSeriesReconcileParams(viewGet(view,'groupName'),[],'scanList',scanList);
  end
end

% Get analysis parameters from averageTSeriesGUI.
nScans = viewGet(view,'nScans');
if ieNotDefined('params')
  % Initialize analysis parameters with default values
  params = averageTSeriesGUI('groupName',viewGet(view,'groupName'));
else
  % Reconcile params with current status of group and ensure that it has
  % the required fields. 
  params = averageTSeriesReconcileParams(params.groupName,params);
end

% Abort if params empty
if ieNotDefined('params')
  mrMsgBox('averageTSeries cancelled');
  return
end

% if just getting params then return
if justGetParams,return,end

% Retrieve parameters
scanList = params.scanList;
reverseList = params.reverseList;
shiftList = params.shiftList;
baseScan = params.baseScan;
groupName = params.groupName;
aveGroupName = params.aveGroupName;
description  = params.description;
interpMethod = params.interpMethod;

% Open new view with the base group
viewBase = newView(viewGet(view,'viewType'));
groupNum = viewGet(viewBase,'groupNum',groupName);
if (groupNum == 0)
  mrErrorDlg('averageTSeries: ',groupName,' does not exist.');
end
viewBase = viewSet(viewBase,'currentGroup',groupNum);

% Open new view and set its group to the average group name. Create the
% group if necessary.
viewAverage = newView(viewGet(view,'viewType'));
aveGroupNum = viewGet(viewAverage,'groupNum',aveGroupName);
if isempty(aveGroupNum)
  view = viewSet(view,'newgroup',aveGroupName);
  aveGroupNum = viewGet(viewAverage,'groupNum',aveGroupName);
end
viewAverage = viewSet(viewAverage,'currentGroup',aveGroupNum);

% Check that all scans in scanList have the same nframes, frameperiod,
% scanvoxelsize, scandims. Note that we check the nFrames and framePeriod
% for the first scan in the list. For other parameters it is relevant
% to which scan we are warping to
nFrames = viewGet(viewBase,'nFrames',scanList(1));
framePeriod = viewGet(viewBase,'framePeriod',scanList(1));
voxelSize = viewGet(viewBase,'scanvoxelsize',baseScan);
scanDims = viewGet(viewBase,'scandims',baseScan);
vol2mag = viewGet(viewBase,'scanVol2mag',baseScan);
vol2tal = viewGet(viewBase,'scanVol2tal',baseScan);
for iscan = 1:length(scanList)
  if (viewGet(viewBase,'nFrames',scanList(iscan)) ~= nFrames)
    mrWarnDlg('(averageTSeries) Can not average these scans because they have different numFrames.');
    keyboard
    return
  end
  if (viewGet(viewBase,'framePeriod',scanList(iscan)) ~= framePeriod)
    mrWarnDlg('These scans  have different frame periods.');
  end
  %     if any(viewGet(viewBase,'scanvoxelsize',scanList(iscan)) ~= voxelSize)
  %         mrErrorDlg('Can not average these scans because they have different voxel sizes.');
  %     end
  if any(viewGet(viewBase,'scandims',scanList(iscan)) ~= scanDims)
    disp(sprintf('(averageTSeries) Scan %i has a different size from the base scan',scanList(iscan)));
  end
end

% Compute output volume
aveTSeries = zeros([scanDims(1) scanDims(2) scanDims(3) nFrames]);
waitHandle = mrWaitBar(0,'Computing average tSeries.  Please wait...');
nSlices = viewGet(view,'nSlices', baseScan);

for iscan = 1:length(scanList)
  scanNum = scanList(iscan);
  reverse = reverseList(iscan);
  shift = shiftList(iscan);

  tseries = loadTSeries(viewBase,scanNum);
  
  % Dump junk frames
  junkFrames = viewGet(viewBase,'junkframes',scanNum);
  tseries = tseries(:,:,:,junkFrames+1:junkFrames+nFrames);

  
  % Time shift
  tseries = circshift(tseries, [0 0 0 shift]);
  
  % Time reverse
  if reverse
    tseries = flipdim(tseries,4);
  end
  
  % Compute transform
  scan2scan = viewGet(view,'scan2scan',baseScan,groupNum,scanNum,groupNum);

  % only warp if needed
  if ~isequal(scan2scan, eye(4))
    % swapXY seems to be needed here, presumably becuase of the way that 
    % warpAffine3 works.
    swapXY = [0 1 0 0;1 0 0 0;0 0 1 0; 0 0 0 1];

    % compute transform
    M = swapXY * scan2scan * swapXY;
    
    % display transformation
    disp(sprintf('Transforming %s:%i to match %s:%i with transformation: ',groupName,scanNum,groupName,baseScan));
    for rownum = 1:4
      disp(sprintf('[%0.2f %0.2f %0.2f %0.2f]',M(rownum,1),M(rownum,2),M(rownum,3),M(rownum,4)));
    end
    % Warp the frames
    waitHandle = mrWaitBar(0,['Warping image volumes for scan ', num2str(iscan)]);
    for frame = 1:nFrames
      mrWaitBar(frame/nFrames,waitHandle);
      tseriesWarped(:,:,:,frame) = warpAffine3(tseries(:,:,:,frame),M,NaN,0,interpMethod,scanDims);
    end  
    % write tseriesWarped back into tseries. This is needed because tseries
    % and tseriesWarped might be different sizes 
    tseries = tseriesWarped;
    clear tseriesWarped;
    mrCloseDlg(waitHandle);
  end
  
  % Add 'em up
  
  for iSlice = 1:nSlices
    tmp = cat(5, aveTSeries(:,:,iSlice,:), tseries(:,:,iSlice,:));
    aveTSeries(:,:,iSlice,:) = nansum(tmp,5);
  end
  % Update waitbar
  %mrWaitBar( iscan/length(scanList) + (iSlice/nSlices * 1/length(scanList)) ,waitHandle);
  mrWaitBar( iscan/length(scanList), waitHandle);

  % remember origianl file/group	
  scanParams.originalFileName{iscan} = viewGet(viewBase,'tSeriesFile',scanNum);
  scanParams.originalGroupName{iscan} = viewGet(viewBase,'groupName',viewGet(viewBase,'curGroup'));
  
end

% Divide by number of scans in scanList
aveTSeries = aveTSeries / length(scanList);
mrCloseDlg(waitHandle);

% Save aveTSeries (using header of 1st scan on scanList as the template for
% the nifti header), and add it as a new scan.
scanParams.fileName = [];
scanParams.junkFrames = 0;
scanParams.totalJunkedFrames = junkFrames;
scanParams.nFrames = nFrames;
scanParams.description = description;
scanParams.vol2mag = vol2mag;
scanParams.vol2tal = vol2tal;
hdr = mlrImageReadNiftiHeader(viewGet(view,'tseriesPath',baseScan));
[viewAverage,tseriesFileName] = saveNewTSeries(viewAverage,aveTSeries,scanParams,hdr);

% Save evalstring for recomputing and params
evalstr = ['view = newView(','''','Volume','''','); view = averageTSeries(view,params);'];
[pathstr,filename] = fileparts(tseriesFileName);
tseriesdir = viewGet(viewAverage,'tseriesdir');
save(fullfile(tseriesdir,filename),'evalstr','params','tseriesFileName');

% Delete temporary viewBase and viewAverage
deleteView(viewBase);
deleteView(viewAverage);

return; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% *** Testing/debugging *** %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params.scanList = 1;
params.shiftList = 2;
params.reverseList = 0;
params.baseScan = 1;
params.groupName = 'Raw';
params.aveGroupName = 'Averages';
params.fileName = [];
params.description = 'test';

MLR.views{1} = averageTSeries(MLR.views{1},params);

MLR.views{1} = viewSet(MLR.views{1},'currentGroup',1);
ts1 = loadTSeries(MLR.views{1},1,6);
ts1 = reshapeTSeries(ts1);
MLR.views{1} = viewSet(MLR.views{1},'currentGroup',2);
ts2 = loadTSeries(MLR.views{1},1,6);
ts2 = reshapeTSeries(ts2);
figure; plot([ts1(9:168,2000) ts2(:,2000)]);



