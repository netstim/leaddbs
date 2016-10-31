% ApplyMotionComptransform.m
%
%        $Id: applyMotionCompTransform.m  2011-12-11 18:57:54Z rosa $
%      usage: thisView = motionCompApplyWarp(thisView,params)
%         by: rosa sanchez and julien besle
%       date: 14/12/2011
%    purpose: applies existing motion compensation parameters estimated on some scans to other scans
%             (useful for applying motion parameters estimated on first echo tp second echo in a dual sequence)
%

function [thisView,params] = applyMotionCompTransform(thisView,params,varargin)

% other arguments
eval(evalargs(varargin));
if ieNotDefined('justGetParams'),justGetParams = 0;end

groupNames=viewGet(thisView,'groupNames');
currentGroupName = viewGet(thisView,'groupName',viewGet(thisView,'currentGroup'));

if ieNotDefined('params')
  params=struct;
end
if fieldIsNotDefined(params,'motionCorrectedGroupName')
  params.motionCorrectedGroupName=groupNames{1};
end
if fieldIsNotDefined(params,'targetsGroupName')
  params.targetsGroupName=currentGroupName;
end
if fieldIsNotDefined(params,'destinationGroupName')
  params.destinationGroupName=groupNames{1};
end
if fieldIsNotDefined(params,'motionCorrectedScanList')
  [~,motionCorrectedGroupNum] = ismember(params.motionCorrectedGroupName,groupNames);
  params.motionCorrectedScanList=1:viewGet(thisView,'nScans',motionCorrectedGroupNum);
end
if fieldIsNotDefined(params,'targetScanList')
  [~,targetsGroupNum] = ismember(params.targetsGroupName,groupNames);
  params.targetScanList=1:viewGet(thisView,'nScans',targetsGroupNum);
end

keepAsking=true;
motionCorrectedGroupNames=groupNames;
targetsGroupNames=groupNames;
destinationGroupNames=groupNames;
while keepAsking
  motionCorrectedGroupNames = putOnTopOfList(params.motionCorrectedGroupName,motionCorrectedGroupNames);
  targetsGroupNames = putOnTopOfList(params.targetsGroupName,targetsGroupNames);
  destinationGroupNames = putOnTopOfList(params.destinationGroupName,destinationGroupNames);
  paramsInfo = {...
    {'motionCorrectedGroupName',motionCorrectedGroupNames,'type=popupmenu','Name of the group of the motion-compensated scans'},...
    {'targetsGroupName',targetsGroupNames,'type=popupmenu','Name of the group of the scans to motion correct (Note: corresponding scans must be identically ordered in the motionCorrected and target groups)'},...
    {'destinationGroupName',destinationGroupNames,'type=popupmenu','Name of the group to put the new motion-compensated scans into'},...
    };
  tempParams = mrParamsDialog(paramsInfo,'applyMotionCompTransform parameters');
  if isempty(tempParams)
    mrMsgBox('motion compensation cancelled');
    return
  else
    [~,motionCorrectedGroupNum] = ismember(tempParams.motionCorrectedGroupName,groupNames);
    if ~strcmp(params.motionCorrectedGroupName,tempParams.motionCorrectedGroupName)
      params.motionCorrectedScanList=1:viewGet(thisView,'nScans',motionCorrectedGroupNum);
    end
    while keepAsking
      tempMotionCompScanList = selectInList(thisView,'scans','Select Motion-compensated scans',params.motionCorrectedScanList,motionCorrectedGroupNum);
      if isempty(tempMotionCompScanList)
        if size(tempMotionCompScanList,2)==1
          mrMsgBox('motion compensation cancelled');
          return;
        else
          break;
        end
      else
        params.motionCorrectedScanList=tempMotionCompScanList;
        [~,targetsGroupNum] = ismember(tempParams.targetsGroupName,groupNames);
        if ~strcmp(params.targetsGroupName,tempParams.targetsGroupName)
          params.targetScanList = 1:viewGet(thisView,'nScans',targetsGroupNum);
        end
        while keepAsking
          tempTargetScanList = selectInList(thisView,'scans','Select target scans',params.targetScanList,targetsGroupNum);
          if isempty(tempTargetScanList)
            if size(tempTargetScanList,2)==1
              mrMsgBox('motion compensation cancelled');
              return;
            else
              break;
            end
          else
            params.targetScanList = tempTargetScanList;
            if length(params.targetScanList)~=length(params.motionCorrectedScanList)
              mrWarnDlg('(applyMotionCompTransform) the number of motion-compensated and target scans must be identical')
            else
              keepAsking=false;
            end
          end
        end
      end
    end
  end
  params = mrParamsCopyFields(tempParams, params);
end

% if just getting params then return
if justGetParams,return,end

% Open new view and set its group to the motion comp group name. Create the
% group if necessary.
viewDestination = newView;
destinationGroupNum = viewGet(viewDestination,'groupNum',params.destinationGroupName);
if isempty(destinationGroupNum)
  thisView = viewSet(thisView,'newgroup',params.destinationGroupName);
  destinationGroupNum = viewGet(viewDestination,'groupNum',params.destinationGroupName);
end
viewDestination = viewSet(viewDestination,'currentGroup',destinationGroupNum);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through scans, warp according to motion estimates and save new tseries%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for s = 1:length(params.targetScanList)
  targetScanNum = params.targetScanList(s);
  transforms = viewGet(thisView,'transforms',params.motionCorrectedScanList(s),motionCorrectedGroupNum);
  totalFrames = viewGet(thisView,'totalFrames',targetScanNum,targetsGroupNum); 
  
  if isempty(transforms)
    mrWarnDlg(['(applyMotionCompTransform) Scan ' num2str(params.motionCorrectedScanList(s)) ' in group ' params.motionCorrectedGroupName ' does not contain motion parameter estimates, skipping...']);
  elseif length(transforms) ~= totalFrames
    mrWarnDlg('(applyMotionCompTransform) Number of frames in motion-corrected and target scans do not match, skipping...');
  elseif any(any(abs(eye(4)- viewGet(thisView, 'scan2scan', params.motionCorrectedScanList(s), motionCorrectedGroupNum, targetScanNum, targetsGroupNum))>1e-6))
    mrWarnDlg('applyMotionCompTransform not implemented for scans with different sforms, skipping...');
  else
    motionParams = getfield(viewGet(thisView,'params',params.motionCorrectedScanList(s),motionCorrectedGroupNum),'params');

    % get the base qform and sform.  the warpedTSeries will inherit these
    baseQform = viewGet(thisView, 'scanqform', motionParams.baseScan, params.targetsGroupName);
    baseSform = viewGet(thisView, 'scansform', motionParams.baseScan, params.targetsGroupName);
    
    % Load tseries
    tseries = loadTSeries(thisView,targetScanNum,'all',[],[],[],[],targetsGroupNum);
    if motionParams.sliceTimeCorrection
      sliceTimes = viewGet(thisView,'sliceTimes',targetScanNum,targetsGroupNum);
    else
      sliceTimes = [];
    end

    % Initialize the warped time series to zeros. Need to re-initialize this
    % for each scan because number of frames can differ. 
    warpedTseries = zeros(size(tseries));

    %this is taken from motionCompPreProcessing 
    % Get slice times and replicate the last frame of tseries for slice time
    % correction
    if motionParams.sliceTimeCorrection
    %  correctedTseries(:,:,:,end+1) = correctedTseries(:,:,:,end); %JB I commented this because the last frame should be replicated only after preprocessing (if it is at all replicated, see end of function) 
      switch motionParams.sliceTimeString
        case 'end of TR'
          newSliceTimes = sliceTimes - 1;  %JB: by removing 1TR to all the acquisition times, 
                                        %the last acquired frame becomes the one with the time closest to 0, hence the reference frame
                                        %(and the first-acquired slices will have to be interpolated closer to their respective next frame)
        case 'middle of TR'             
          newSliceTimes = sliceTimes - 0.5; %JB: same with .5 TR
        case 'beginning of TR'         
          %JB: the first acquired frame is the reference time
          % so the slice acquisition times do not change (0 for the first acquired frame and close to 1 for the last-acquired frame,
          % which will have to be interpolated closer to the previous frame)
          newSliceTimes= sliceTimes;
        otherwise
          mrErrorDlg('Invalid slice times');
      end
    else
      newSliceTimes = [];
    end

    % warp the images according to the motion estimates
    waitHandle = mrWaitBar(0,['(applyMotionCompTransform) Warping image volumes for scan ',num2str(targetScanNum),'.  Please wait...']);
    for frame = 1:totalFrames
      mrWaitBar(frame/totalFrames,waitHandle)
      if motionParams.sliceTimeCorrection
        warpedTseries(:,:,:,frame) = warpAffineInterp3(tseries,frame,transforms{frame},newSliceTimes,NaN,motionParams.interpMethod);
      else
        warpedTseries(:,:,:,frame) = warpAffine3(tseries(:,:,:,frame),transforms{frame},NaN,0,motionParams.interpMethod);
      end
    end
    mrCloseDlg(waitHandle);
    clear tseries;

    % Save tseries with modified nifti header
    scanParams = viewGet(thisView,'scanParams',targetScanNum,targetsGroupNum);
    scanParams.description = ['Full Motion compensation of ',params.targetsGroupName,' scan ',num2str(targetScanNum),': ',scanParams.description];
    scanParams.fileName = [];
    scanParams.originalFileName{1} = viewGet(thisView,'tseriesfile',targetScanNum,targetsGroupNum);
    scanParams.originalGroupName{1} = viewGet(thisView,'groupName',targetsGroupNum);

    % set the qform and sform to that of the baseScan
    scanParams.niftiHdr = cbiSetNiftiQform(scanParams.niftiHdr, baseQform);
    scanParams.niftiHdr = cbiSetNiftiSform(scanParams.niftiHdr, baseSform);

    [viewDestination,tseriesFileName] = saveNewTSeries(viewDestination,warpedTseries,scanParams,scanParams.niftiHdr);

    % Save evalstring for recomputing and params
    evalstr = ['thisView = newView(','''','Volume','''','); thisView = applyMotionCompTransform(thisView,params);'];
    [pathstr,filename] = fileparts(tseriesFileName);
    tseriesdir = viewGet(viewDestination,'tseriesdir');
    save(fullfile(tseriesdir,filename),'evalstr','params','transforms','tseriesFileName');

    % clear temporary tseries
    clear warpedTseries;
  end
end

% Delete temporary views
deleteView(viewDestination);

return

