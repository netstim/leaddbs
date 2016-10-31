% fslFugueTSeries.m
%
%       $Id: fslFugueTSeries.m 1833 2010-11-13 18:37:14Z julien $	
%      usage: fslFugueTSeries(view,params,varargin)
%         by: julien besle
%       date: 11/05/2012
%    purpose: Undistort EPI time-series using FSL FUGUE
%
% Examples:
%
% params = fslFugueTSeries('groupName','Raw');
% view = newView;
% view = fslFugueTSeries(view,params);
%
% view = fslFugueTSeries(view);
%
% To just get parameters
% [view params] = fslFugueTSeries(view,[],'justGetParams=1','defaultParams=1','scanList=[1 8]');


function [view params] = fslFugueTSeries(view,params,varargin)

if strcmp(mrGetPref('fslPath'),'FSL not installed')
  mrWarnDlg('(fslApplyWarpOverlays) No path was provided for FSL. Please set MR preference ''fslPath'' by running mrSetPref(''fslPath'',''yourpath'')')
  return
end

% other arguments
eval(evalargs(varargin));
if ieNotDefined('justGetParams'),justGetParams = 0;end
if ieNotDefined('defaultParams'),defaultParams = 0;end
if ieNotDefined('scanList'),scanList = [];end

nScans = viewGet(view,'nScans');
if (nScans == 0)
  mrWarnDlg('(fslFugueTSeries) No scans in group');
  return
end

% Get analysis parameters from fslFugueTSeriesGUI.
if ieNotDefined('params')
  % Initialize analysis parameters with default values
  params = fslFugueTSeriesGUI('groupName',viewGet(view,'groupName'),'defaultParams',defaultParams,'scanList',scanList,'v',view);
else
  % Reconcile params with current status of group and ensure that it has
  % the required fields.
  params = fslFugueTSeriesReconcileParams(params.groupName,params);
end

% Abort if params empty
if ieNotDefined('params')
  mrMsgBox('FSL FUGUE undistortion cancelled');
  return
end

% if just getting params then return
if justGetParams,return,end

% Retrieve parameters
fieldMapFile = params.fieldMapFile;
dwellTime = params.dwellTime;
unwarpDirection = params.unwarpDirection;
targetScans = params.targetScans;
undistortedGroupName = params.undistortedGroupName;
groupName = params.groupName;
descriptions  = params.descriptions;
useMask=params.useMask;
multiplyBy2pi = params.multiplyBy2pi;

% Open new view with the base group
baseView = newView;
groupNum = viewGet(baseView,'groupNum',groupName);
if (groupNum == 0)
  mrErrorDlg('fslFugueTSeries: ',groupName,' does not exist.');
end
baseView = viewSet(baseView,'currentGroup',groupNum);


% Ignore scans if data size is different from that for base scan.
% *** get rid of this by using warpAffine3 below to resample each scan to base scan size
hdr=mlrImageReadNiftiHeader(params.fieldMapFile);
for scanNum = targetScans
  datasize = viewGet(baseView,'datasize',scanNum);
  if any(datasize ~= hdr.dim(2:4)')
    mrWarnDlg(['(fslFugueTSeries) Ignoring scan ',num2str(scanNum),'. The scan dimensions must match the field map dimensions']);
    targetScans = targetScans(targetScans ~= scanNum);
  end
end

% Open new view and set its group to the motion comp group name. Create the
% group if necessary.
fslFugueView = newView;
undistortedGroupNum = viewGet(fslFugueView,'groupNum',undistortedGroupName);
if isempty(undistortedGroupNum)
  view = viewSet(view,'newgroup',undistortedGroupName);
  undistortedGroupNum = viewGet(fslFugueView,'groupNum',undistortedGroupName);
end
fslFugueView = viewSet(fslFugueView,'currentGroup',undistortedGroupNum);

tempFileName = 'temp.hdr';
tempFieldMapFileName = 'temp2.hdr';
B0correctedTseriesdir = viewGet(fslFugueView,'tseriesdir');
tseriesdir = viewGet(baseView,'tseriesdir');

if useMask
  [mask view] = motionCompGetMask(view,params,scanNum,groupNum);
  cbiWriteNifti(tempFieldMapFileName, mask, hdr);
  %mask the field map with the ROI mask using fslmaths
  command =  sprintf('fslmaths %s -mas %s %s',...
  fieldMapFile, tempFieldMapFileName, tempFieldMapFileName);
  executeCommand(command, 'fslmaths -mas')
  fieldMapFile=tempFieldMapFileName;
end

if multiplyBy2pi
  %multiply the field map by 2*pi using fslmaths
  command =  sprintf('fslmaths %s -mul 6.283185307179586 %s',...
  fieldMapFile,tempFieldMapFileName);
  executeCommand(command, 'fslmaths -mul')
  fieldMapFile=tempFieldMapFileName;
end

for iScan = 1:length(targetScans)
  scanNum = targetScans(iScan);
  scanFileName =viewGet(baseView,'tseriesfile',scanNum);
  % get the base qform and sform.  the warpedTSeries will inheret these
  baseQform = viewGet(baseView, 'scanqform', scanNum, groupName);
  baseSform = viewGet(baseView, 'scansform', scanNum, groupName);
  
  
  command =  sprintf('fugue -i %s --loadfmap=%s --dwell=%f -u %s --unwarpdir=%s',...
    fullfile(tseriesdir,scanFileName),fieldMapFile, dwellTime, tempFileName, unwarpDirection);
  executeCommand(command, 'fugue')
  
  % Save tseries with modified nifti header
  scanParams = viewGet(baseView,'scanParams',scanNum);
  scanParams.description = descriptions{iScan};
  scanParams.fileName = [];
  scanParams.originalFileName{1} = scanFileName;
  scanParams.originalGroupName{1} = viewGet(baseView,'groupName');

  % set the qform and sform to that of the dwellTime
  scanParams.niftiHdr = cbiSetNiftiQform(scanParams.niftiHdr, baseQform);
  scanParams.niftiHdr = cbiSetNiftiSform(scanParams.niftiHdr, baseSform);
  
  [fslFugueView,tseriesFileName] = saveNewTSeries(fslFugueView,tempFileName,scanParams,scanParams.niftiHdr);

  %if frame period was in milliseconds in original file, fugue changes it to second without changing the value, fix that
  newFramePeriod=viewGet(fslFugueView,'framePeriod',viewGet(fslFugueView,'nscans'));
  if newFramePeriod~=scanParams.framePeriod
    newHdr=mlrImageReadNiftiHeader(fullfile(B0correctedTseriesdir,tseriesFileName));
    niftiSpaceUnit = rem(newHdr.xyzt_units, 8); 
    niftiTimeUnit = rem(newHdr.xyzt_units-niftiSpaceUnit, 64);
    switch(niftiTimeUnit)
      case 8
        timeUnit = 'sec';
      case 16
        timeUnit = 'msec';
      case 32
        timeUnit = 'microsec';
    end
    fprintf('(fslfugueTSeries) Changing frame period from %f %s to %f sec in file %s\n',newFramePeriod,timeUnit,scanParams.framePeriod,tseriesFileName);
    newHdr.pixdim(5)=scanParams.framePeriod;
    cbiWriteNiftiHeader(newHdr,fullfile(B0correctedTseriesdir,tseriesFileName));
    scanParams.fileName=tseriesFileName;
    viewSet(fslFugueView,'updateScan',scanParams, viewGet(fslFugueView,'nscans'));
  end

  % Save evalstring for recomputing and params
  evalstr = ['view = newView(','''','Volume','''','); view = fslFugueTSeries(view,params);'];
  [~,filename] = fileparts(tseriesFileName);
  save(fullfile(B0correctedTseriesdir,filename),'evalstr','params');

end

disp('(fslfugueTSeries) Done');

% Delete temporary views
deleteView(baseView);
deleteView(fslFugueView);

return

function executeCommand(command, commandName)

disp(command);
try
  [s,w] = unix(command);

  if s ~= 0 % unix error
    disp('UNIX error message:')
    disp(w)
    disp('-------------------')
    return
  end
catch 
  fprintf('(fslApplyWarp) There was a problem running FSL %s command',commandName)
  disp(sprintf('unix error code: %d; %s', s, w))
  return
end
