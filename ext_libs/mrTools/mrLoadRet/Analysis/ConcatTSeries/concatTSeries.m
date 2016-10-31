% concatTSeries.m
%
%        $Id$
%      usage: concatTSeries(view, params)
%         by: justin gardner
%       date: 10/12/06
%    purpose: concatenate time series together.
%
%             to just get a default parameter structure:
% 
%             v = newView;
%             [v params] = concatTSeries(v,[],'justGetParams=1');
%             [v params] = concatTSeries(v,[],'justGetParams=1','defaultParams=1');
%             [v params] = concatTSeries(v,[],'justGetParams=1','defaultParams=1','scanList=[1 2]');
%
%
function [view params] = concatTSeries(view,params,varargin)

% check arguments
if ~any(nargin == [1 2 3 4 5 6 7 8])
  help concatTSeries.m
  return
end

% other arguments
eval(evalargs(varargin));
if ieNotDefined('justGetParams'),justGetParams = 0;end
if ieNotDefined('defaultParams'),defaultParams = 0;end

interpTypes = {'nearest','linear','spline','cubic'};
if find(strcmp(mrGetPref('interpMethod'),interpTypes))
  interpTypes = putOnTopOfList(mrGetPref('interpMethod'),interpTypes);
end

% check to see if any of the scans in this group have a non identity scan2scan
needToWarp = 0;
for i = 2:viewGet(view,'nScans')
  if ~isequal(viewGet(view,'scan2scan',1,[],i),eye(4))
    needToWarp = 1;
  end
end

% check to see if any scans have a tSense that is not one, if it
% is then we will offer to notch filter the data, otherwise hide the option
offerNotchFilter = false;
defaultNotchFilterSetting = false;
if ~isempty(viewGet(view,'groupNum','Raw'))
  for iScan = 1:viewGet(view,'nScans','Raw')
    tSense = viewGet(view,'auxParam','tSense',iScan,'Raw');
    if iscell(tSense),tSense = cell2mat(tSense);end
    if isscalar(tSense) && (tSense > 1)
      offerNotchFilter = true;
      % set default setting to false if the current scan is a concat
      if ~isempty(viewGet(view,'concatInfo'))
        defaultNotchFilterSetting = false;
      end
      break;
    end
  end
end

% description of paramaters (used by mrParamsDialog functions)
paramsInfo = {...
    {'groupName',putOnTopOfList(viewGet(view,'groupName'),viewGet(view,'groupNames')),'Name of group from which to make concatenation'},...
    {'newGroupName','Concatenation','Name of group to which concatenation will be saved. If group does not already exist, it will be created.'},...
    {'description','Concatenation of [x...x]','Description that will be set to have the scannumbers that are selected'},...
    {'filterType',{'Detrend and highpass','Detrend only','None'},'Which filter to use. Highpass filtering will use the cutoff below.'},...
    {'filterCutoff',0.01,'minmax=[0 inf]','Highpass filter cutoff in Hz'},...
    {'percentSignal',1,'type=checkbox','Convert to percent signal change. This is done by simply dividing by the mean, so that you get a timecourse where the mean is 1. (The mean is not subtracted out so that if subsequent analysis tries to divide by mean again it will not affect the time series)'}};
if needToWarp
  paramsInfo{end+1} = {'warp',1,'type=checkbox','Warp images based on alignment. This can be used to concatenate together scans taken on different days. If the scans have the same slice prescription this will not do anything.'};
  paramsInfo{end+1} = {'warpInterpMethod',interpTypes,'Interpolation method for warp','contingent=warp'};
end
paramsInfo{end+1} = {'projectOutMeanVector',0,'type=checkbox','Project out a mean vector defined from one roi out of the data. This is used if you want to remove a global component that is estimated as the mean over some roi from your data. If you select this you will be asked to choose one roi for defining the mean vector of what you want to project out and another roi from which you want to project out.'};
if offerNotchFilter
  paramsInfo{end+1} = {'notchFilterForTSense',defaultNotchFilterSetting,'type=checkbox','This is used to notch out the highest frequency for tSense data'};
end
    
% First get parameters
if ieNotDefined('params')
  % Initialize analysis parameters with default values
  if defaultParams
    params = mrParamsDefault(paramsInfo);
  else
    params = mrParamsDialog(paramsInfo);
  end
  % no params means user hit cancel
  if isempty(params),return,end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % select scans
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  view = viewSet(view, 'groupName', params.groupName);
  if ~ieNotDefined('scanList')
    params.scanList = scanList;
  elseif defaultParams
    params.scanList = 1:viewGet(view,'nScans');
  else
    params.scanList = selectInList(view,'scans','Select Scans',1:viewGet(view,'nScans'));
  end
  if isempty(params.scanList),return,end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % if warp is set, then ask which scan to use as base scan for warp, unless set by default
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if ~needToWarp,params.warp = 0;params.warpInterpMethod = interpTypes{1};end
  if params.warp
    if defaultParams
      params.warpBaseScan = params.scanList(1);
    else
      % create a list of scans to ask user which scan to use as the base scan
      for i = 1:viewGet(view,'nScans')
	scanNames{i} = sprintf('%i:%s (%s)',i,viewGet(view,'description',i),viewGet(view,'tSeriesFile',i));
      end
      warpParams = mrParamsDialog({{'warpBaseScan',scanNames,'The scan that will be used as the base scan to warp all the other scans to'}});
      if isempty(warpParams),return,end
      params.warpBaseScan = find(strcmp(warpParams.warpBaseScan,scanNames));
    end
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Get rois for doing projection
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if params.projectOutMeanVector || defaultParams
    % projectOutMeanVectorParams
    params.projectOutMeanVectorParams = projectOutMeanVectorParams(view,defaultParams);
    if isempty(params.projectOutMeanVectorParams),return,end
  end
  % check the parameters
  params = mrParamsReconcile(params.groupName,params);
else
  % Reconcile params with current status of group and ensure that it has
  % the required fields. 
  params.paramsInfo = paramsInfo;
  params = mrParamsReconcile(params.groupName,params);
end
drawnow;

% Abort if params empty
if ieNotDefined('params'),return,end

% if just getting params then return
if justGetParams,return,end
  
% Open new view with the base group
viewBase = newView;
groupNum = viewGet(viewBase,'groupNum',params.groupName);

if (groupNum == 0)
  mrErrorDlg('concatTSeries: ',groupName,' does not exist.');
end
viewBase = viewSet(viewBase,'currentGroup',groupNum);

% Open new view and set its group to the concat group name. Create the
% group if necessary.
viewConcat = newView;
concatGroupNum = viewGet(viewConcat,'groupNum',params.newGroupName);
if isempty(concatGroupNum)
  view = viewSet(view,'newgroup',params.newGroupName);
  concatGroupNum = viewGet(viewConcat,'groupNum',params.newGroupName);
end
viewConcat = viewSet(viewConcat,'currentGroup',concatGroupNum);

% Check that all scans in scanList have the same 
% scanvoxelsize, scandims

% get the base scan, if we are doing a warp then choose
% the warpBaseScan, otherwise choose the first scan in the scan list
if ~params.warp
  baseScanNum = params.scanList(1);
else
  baseScanNum = params.warpBaseScan;
end

% get info about scan
d.tr = viewGet(viewBase,'framePeriod',baseScanNum);
d.voxelSize = viewGet(viewBase,'scanvoxelsize',baseScanNum);
d.dim = viewGet(viewBase,'scandims',baseScanNum);
d.nFrames = viewGet(viewBase,'nFrames',baseScanNum);
d.dim(4) = d.nFrames;
d.sform = viewGet(viewBase,'scanSform',baseScanNum);
d.vol2mag = viewGet(viewBase,'scanVol2mag',baseScanNum);
d.vol2tal = viewGet(viewBase,'scanVol2tal',baseScanNum);

for iscan = 1:length(params.scanList)
  % check frame period for mismatch
  if (viewGet(viewBase,'framePeriod',params.scanList(iscan)) ~= d.tr)
    mrWarnDlg(sprintf('concatTSeries: These scans have different TR. (%0.4f vs %0.4f)',viewGet(viewBase,'framePeriod',params.scanList(iscan)),d.tr));
  end
  % make sure that the voxel sizes are the same to within roundoff error
  baseVoxelSize = viewGet(viewBase,'scanvoxelsize',params.scanList(iscan));
  roundoff = 100000;
  if ~isequal(round(baseVoxelSize*roundoff)/roundoff,round(d.voxelSize*roundoff)/roundoff)
    disp(sprintf('(concatTSeries) Scans have different voxel sizes %i:[%s]~=[%s]',params.scanList(iscan),num2str(baseVoxelSize),num2str(d.voxelSize)));
  end
  % check the scan dims
  if ~isequal(viewGet(viewBase,'scandims',params.scanList(iscan)),d.dim(1:3))
    if params.warp
      disp('(concatTSeries) Scans have different dimensions.');
    else
      mrErrorDlg('(concatTSeries) Must choose warp in order to concat scans with different sizes.');
    end
  end
  % check the scan sforms
  if ~isequal(viewGet(viewBase,'scanSform',params.scanList(iscan)),d.sform)
    % this is only an issue if warp is  not set
    if ~params.warp
      mrWarnDlg(sprintf('(concatTSeries) Sform for scan %s:%i does not match %s:%i. This means that they have different slice prescriptions. Usually you should select warp in this case so that the different scans are all warped together to look like the base scan. You have not selected warp.',params.groupName,params.scanList(iscan),params.groupName,baseScanNum));
    end
  end
  % check the vol2mag and vol2tal
  if (~isequal(viewGet(viewBase,'scanVol2mag',params.scanList(iscan)),d.vol2mag) || ...
      ~isequal(viewGet(viewBase,'scanVol2tal',params.scanList(iscan)),d.vol2tal))
    % this is only an issue if warp is  not set
    if ~params.warp
      mrWarnDlg(sprintf('(concatTSeries) The scanVol2mag/scanVol2tal for scan %s:%i does not match %s:%i. This means that they have been aligned to different volume anatomies. Usually you should select warp in this case so that the different scans are all warped together to look like the base scan. You have not selected warp.',params.groupName,params.scanList(iscan),params.groupName,baseScanNum));
    end
  end
end
disp(sprintf('(concatTSeries) FramePeriod for scan is: %0.2f',d.tr));

% initialize some things
concatInfo.n = 0;
concatInfo.whichScan = [];
concatInfo.whichVolume = [];

set(viewGet(view,'figNum'),'Pointer','watch');drawnow;
tic
% Compute output volume
waitHandle = mrWaitBar(0,'(concatTSeries) Concatenating tSeries.  Please wait...');
view = viewSet(view,'curGroup',groupNum);
for iscan = 1:length(params.scanList)
  scanNum = params.scanList(iscan);

  viewBase = viewSet(viewBase,'curScan',params.scanList(iscan));
  view = viewSet(view,'curScan',scanNum);
  
  % Load it
  mrDisp(sprintf('\n(concatTSeries) Loading scan %i from %s\n',scanNum,viewGet(viewBase,'groupName')));
  tSeries = loadTSeries(viewBase,scanNum,'all');
	
  % Dump junk frames
  junkFrames = viewGet(viewBase,'junkframes',scanNum);
  d.nFrames = viewGet(viewBase,'nFrames',scanNum);
  d.dim(4) = d.nFrames;
  tSeries = tSeries(:,:,:,junkFrames+1:junkFrames+d.nFrames);

  % get the total junked frames. This is the number of frames
  % we have junked here, plus whatever has been junked in previous ones
  thisTotalJunkedFrames = viewGet(viewBase,'totalJunkedFrames',scanNum);
  if length(thisTotalJunkedFrames > 1)
    % give warning if there are non-zero total junked frames
    if sum(thisTotalJunkedFrames) ~= 0
      disp(sprintf('(concatTSeries) This scan has multiple total junked frames [%s] -- i.e. it looks like an average or a concat. Using a junk frame count of %i',num2str(thisTotalJunkedFrames),median(thisTotalJunkedFrames)));
    end
    thisTotalJunkedFrames = median(thisTotalJunkedFrames);
  end    
  totalJunkedFrames(iscan) = junkFrames+thisTotalJunkedFrames;
  
  % Compute transform
  d.data = [];
  if params.warp
    % get the scan2scan xform. Note that we get the scan coordinates form the base scan to the
    % current scan. When we do interpolation we want to know for each destination location (warp base)
    % where the voxels live in the source (current scan). So we need a transform that takes each
    % destination location (warp base scan coordinate) to locations in the current scan.
    scan2scan = viewGet(viewBase,'scan2scan',params.warpBaseScan,groupNum,scanNum,groupNum);
    
    if ~isequal(scan2scan,eye(4))
      % swapXY seems to be needed here, presumably becuase of the way that 
      % warpAffine3 works.
      swapXY = [0 1 0 0;1 0 0 0;0 0 1 0; 0 0 0 1];

      % compute transformation matrix
      M = swapXY * scan2scan * swapXY;

      % display transformation
      disp(sprintf('Transforming %s:%i to match %s:%i with transformation: ',params.groupName,scanNum,params.groupName,params.warpBaseScan));
      for rownum = 1:4
	disp(sprintf('[%0.2f %0.2f %0.2f %0.2f]',M(rownum,1),M(rownum,2),M(rownum,3),M(rownum,4)));
      end
      % Warp the frames
      for frame = 1:d.nFrames
	d.data(:,:,:,frame) = warpAffine3(tSeries(:,:,:,frame),M,NaN,0,params.warpInterpMethod,d.dim);
      end 
    else
      d.data = tSeries; % if didn't warp, just set data
    end
  else
    d.data = tSeries; % if didn't warp, just set data
  end

  % set up notch filter for tSense - this is intended to remove the highest
  % temporal frequency component from the data that are introduced by
  % the tSense reconstruction
  if ~isfield(params,'notchFilterForTSense'),params.notchFilterForTSense = 0;end
  d.notchFilterForTSense = 0;
  if params.notchFilterForTSense
    if ~strcmp(lower(params.filterType),lower('Detrend and highpass'))
      mrWarnDlg(sprintf('(concatTSeries) !!! Notch filter for tSense is only implemented with the highpass filter for now. !!! Not applying any notch filter'));
    else
      % get tSense
      tSense = viewGet(view,'auxParam','tSense');
      if iscell(tSense),tSense = cell2mat(tSense);end
      if isscalar(tSense)
	if any(tSense == [2 4])
	  d.notchFilterForTSense = tSense;
	elseif tSense ~= 1
	  mrWarnDlg(sprintf('(concatTSeries) !!! Notch filter for tSense has only be implemented for acceleration factors of 2 (tSense=%i) !!! Not running any notch filter',tSense));
	end
      else
	disp(sprintf('(concatTSeries) Notch filter not being run on concat of a concat'));
      end
    end
  end
  
  % do other  processing here
  if strcmp(lower(params.filterType),lower('Detrend and highpass'))
    d = eventRelatedHighpass(d,params.filterCutoff);
  elseif strcmp(lower(params.filterType),lower('Detrend only'))
    d = detrendTSeries(d);
  elseif ~strcmp(lower(params.filterType),lower('None'))
    disp(sprintf('(concatTSeries) Unknown filterType; %s, not applying any filtering',params.filterType));
  end
  
  if params.projectOutMeanVector
    % project out the mena vector calculated from the selected roi
    [projection d.data] = projectOutMeanVector(view,params.projectOutMeanVectorParams,[],d.data);
    % then keep some of the information about the projection
    d.projection.sourceName = projection.sourceName;
    d.projection.sourceMeanVector = projection.sourceMeanVector;
    % and the information about the projection magnitude
    d.projection.r = projection.r;
    d.projection.reconProjectionMagnitude = projection.reconProjectionMagnitude;
    d.projection.linearCoords = projection.linearCoords;
    clear projection;
  end
    
  % convert to percent signal change
  warning off                           % to avoid divide by zero warnings...
  if params.percentSignal
    d.mean = mean(d.data,4);
    % for means that are zero, divide by nan
    d.mean(d.mean==0) = nan;

    disppercent(-inf, '(concatTSeries) Converting to percent signal change');
    for i = 1:d.dim(4)
      d.data(:,:,:,i) = (d.data(:,:,:,i)./d.mean);
      if params.percentSignal == 2           % scale it to mean of 1,000
          params.scaleFactor = 10000;
          d.data(:,:,:,i) = d.data(:,:,:,i) * params.scaleFactor;
      end
      disppercent(i/d.dim(4));
    end
    disppercent(inf);
  end
  warning on

  % get the path and filename
  [path,filename,ext] = fileparts(viewGet(viewBase,'tseriesPath',scanNum));
  baseGroupName = viewGet(viewBase,'groupName');

  % Save TSeries (using header of 1st scan on scanList as the template for
  % the nifti header), and add it as a new scan.
  if iscan == 1
    scanParams.fileName = [];
    scanParams.junkFrames = 0;
    scanParams.nFrames = d.nFrames;
    scanParams.description = params.description;
    scanParams.originalFileName{1} = filename;
    scanParams.originalGroupName{1} = baseGroupName;
    scanParams.totalJunkedFrames = totalJunkedFrames;
    scanParams.vol2mag = d.vol2mag;
    scanParams.vol2tal = d.vol2tal;

    hdr = mlrImageReadNiftiHeader(viewGet(viewBase,'tseriesPath',baseScanNum));
    % data *MUST* be written out as float32 b/c of the small values-epm
    hdr.datatype = 16;
    % if we are warping, then we need to change the sform to the
    % sform of the scan we warped to
    if params.warp
      hdr.sform44 = viewGet(viewBase,'scanSform',params.warpBaseScan,groupNum);
      hdr.sform_code = viewGet(viewBase,'scanSformCode',params.warpBaseScan,groupNum);
    end

    [viewConcat,tseriesFileName] = saveNewTSeries(viewConcat,d.data,scanParams,hdr);
    % get new scan number
    saveScanNum = viewGet(viewConcat,'nScans');
    
    % now load up channels
    stimfile = viewGet(viewBase,'stimFile',scanNum);
  % for subsequent scans, we are going to append
  else
    oldScanParams = viewGet(viewConcat,'scanParams',saveScanNum);
    oldScanParams.originalFileName{end+1} = filename;
    oldScanParams.originalGroupName{end+1} = baseGroupName;
    oldScanParams.totalJunkedFrames = totalJunkedFrames;
    viewConcat = saveTSeries(viewConcat,d.data,saveScanNum,oldScanParams,[],1);

  end   
  % remember some info about where data comes from
  concatInfo.n = concatInfo.n+1;
  concatInfo.nifti(concatInfo.n) = hdr;

  % save the original path and filename
  concatInfo.filename{concatInfo.n} = filename;
  concatInfo.path{concatInfo.n} = path;

  % save which scan and volume every frame is from
  concatInfo.whichScan = [concatInfo.whichScan iscan*ones(1,d.dim(4))];
  concatInfo.whichVolume = [concatInfo.whichVolume 1:d.dim(4)];

  % save the junk frames
  concatInfo.junkFrames(concatInfo.n) = junkFrames;
  
  % keep the highpass filter used
  if isfield(d,'hipassfilter')
    concatInfo.hipassfilter{concatInfo.n} = d.hipassfilter;
  end
  concatInfo.filterType = params.filterType;

  % keep the projection info if that was used
  if isfield(d,'projection')
    projectionConcat{concatInfo.n} = d.projection;
    concatInfo.projection{concatInfo.n}.sourceName = d.projection.sourceName;
    concatInfo.projection{concatInfo.n}.sourceMeanVector = d.projection.sourceMeanVector;
  end
    
  % update the runtransition field, which keeps the beginning
  % and end volume of each run
  totalVolumes = length(concatInfo.whichVolume);
  concatInfo.runTransition(iscan,1) = totalVolumes-d.dim(4)+1;
  concatInfo.runTransition(iscan,2) = totalVolumes;

  % Update waitbar
  mrWaitBar(iscan/length(params.scanList),waitHandle);
end
mrCloseDlg(waitHandle);
toc;

% if we have done projection, then save an overlay
if params.projectOutMeanVector
  % make overlay for each scan in concatenation
  meanOverlay = zeros(viewGet(viewConcat,'scanDims',saveScanNum));
  for i = 1:length(projectionConcat)
    % create an overlay
    overlay{i} = zeros(viewGet(viewConcat,'scanDims',saveScanNum));
    overlay{i}(projectionConcat{i}.linearCoords) = projectionConcat{i}.r;
    % make a mean overlay as well
    meanOverlay = meanOverlay+overlay{i}/length(projectionConcat);
    % make a name for the overlay
    overlayNames{i} = sprintf('r%i',i);
    % remove the field since it now stored in the overlay
    projectionConcat{i} = rmfield(projectionConcat{i},'r');
  end
  if length(projectionConcat) > 1
    % add the mean overlay
    overlay{end+1} = meanOverlay;
    overlayNames{end+1} = sprintf('mean_r');
  else
    overlayNames{1} = 'r';
  end
  % and save an analysis with the overlays
  mrDispOverlay(overlay,saveScanNum,viewGet(viewConcat, 'curGroup'),viewConcat,'overlayNames',overlayNames,'saveName=projectionAnal','analName=projectionAnal','range',[-1 1],'clip',[0.1 -0.1],'cmap',[flipud(fliplr(hot(128)));hot(128)],'colormapType','normal','d',projectionConcat,'colormapType=setRangeToMaxAroundZero','interrogator=projectOutMeanVectorPlot');
  % check to see if projection analysis is loaded in concatenation group.
  % if it is unload it and reload it
  curMLRGroup = viewGet(view,'curGroup');
  view = viewSet(view,'curGroup',concatGroupNum)
  analysisNames = viewGet(view,'analysisNames');
  if any(strcmp('projectionAnal',analysisNames))
    disp(sprintf('(concatTSeries) Projection anal is loaded in projection group. Reloading with current analysis'));
    projectionAnalysisNum = find(strcmp('projectionAnal',analysisNames));
    view = viewSet(view,'deleteAnalysis',projectionAnalysisNum);
    view = loadAnalysis(view,fullfile('projectionAnal','projectionAnal.mat'));
  end
  % reset back to original group
  view = viewSet(view,'curGroup',curMLRGroup);
end

% Save evalstring for recomputing and params
evalstr = ['view = newView; view = concatTSeries(view,params);'];
tseriesdir = viewGet(viewConcat,'tseriesdir');
[pathstr,filename,ext] = fileparts(fullfile(tseriesdir,tseriesFileName));
save(fullfile(pathstr,filename),'evalstr','params','concatInfo');

% Delete temporary viewBase and viewConcat
deleteView(viewBase);
deleteView(viewConcat);

set(viewGet(view,'figNum'),'Pointer','arrow');drawnow

return; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% concat channels together
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function concatInfo = concatTraces(concatInfo, newTraces)

% find last acquistition pulse in data set
lastacq = last(find(concatInfo.traces(1,:)));
% find spacing of acq pulses
acqspace = median(diff(find(concatInfo.traces(1,:))));
% find out where the new acq pulses start
firstacq = first(find(newTraces(1,:)));
% find loc of new acq pulses
newacqpos = firstacq:length(newTraces(1,:));
% find length of new acq pulses
newacqlen = length(newacqpos);
% figure out where to put new acq pulses
putnewacqpos = (lastacq+acqspace):(lastacq+acqspace+newacqlen-1);
% put the stim traces there too
concatInfo.traces(:,putnewacqpos) = newTraces(:,newacqpos);


%%%%%%%%%%%%%%%%%%%%%%%%%
%%   detrend tSeries   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function d = detrendTSeries(d)

disppercent(-inf,sprintf('(concatTSeries:detrendTSeries) Detrending data'));
d.data = reshape(eventRelatedDetrend(reshape(d.data,prod(d.dim(1:3)),d.dim(4))')',d.dim(1),d.dim(2),d.dim(3),d.dim(4));
disppercent(inf);





 
