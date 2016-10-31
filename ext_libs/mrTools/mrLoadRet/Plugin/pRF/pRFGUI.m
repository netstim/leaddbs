% pRFGUI.m
%
%        $Id:$ 
%      usage: pRFGUI()
%         by: justin gardner
%       date: 11/20/11
%    purpose: GUI for getting params for pRF
%
function params = pRFGUI(varargin)

% get the arguments
params=[];groupNum=[];defaultParams=[];scanList = [];v = [];pRFFitParamsOnly=[];
getArgs(varargin,{'params=[]','groupNum=[]','defaultParams=0','scanList=[]','v=[]','pRFFitParamsOnly=0'});

% if called with params, then just display
if ~isempty(params)
  retval = dispParams(params);
  if isempty(retval),params = [];end
  return
end

% get a view
deleteViewOnExit = false;
if isempty(v),v = newView;deleteViewOnExit = true;end

% get the group names put on top passed in group if set
groupNames = putOnTopOfList(viewGet(v,'groupName',groupNum),viewGet(v,'groupNames'));

% get possible restrictions of the analysis (i.e. restrict only to volumes
% in base coordinates, rois, etc.)
restrict = {'None'};putOnTop='';nBases = 0;
curBase = viewGet(v,'curbase');
for iBase = 1:viewGet(v,'numbase')
  b = viewGet(v,'base',iBase);
  if b.type > 0
    % add any surface or flat map to base list
    restrict{end+1} = sprintf('Base: %s',b.name);
    if isequal(iBase,curBase)
      putOnTop{1} = restrict{end};
    end
    % update number of bases we have found
    nBases = nBases+1;
  end
end
% if we have found more than 1 base then add the possibility of
% running on all bases
if nBases >= 2
  restrict{end+1} = sprintf('Base: ALL');
end
% get rois
roiNames = viewGet(v,'roiNames');
curROI = viewGet(v,'curROI');
for iRoi = 1:length(roiNames)
  restrict{end+1} = sprintf('ROI: %s',roiNames{iRoi});
  if iRoi == curROI
    putOnTop{end+1} = restrict{end};
  end
end
if ~isempty(putOnTop)
  for i = 1:length(putOnTop)
    restrict = putOnTopOfList(putOnTop{i},restrict);
  end
end

% check if we have an old pRF analysis loaded which
% the user might want to continue running (i.e. add voxels to)
analysisNames = viewGet(v,'analysisNames');
pRFAnalyses = {};
for i = 1:length(analysisNames)
  if isequal(viewGet(v,'analysisType',i),'pRFAnal')
    pRFAnalyses{end+1} = analysisNames{i};
  end
  % put current analysis on top of list
  pRFAnalyses = putOnTopOfList(viewGet(v,'analysisName'),pRFAnalyses);
end

% set the parameter string
paramsInfo = {};
if ~pRFFitParamsOnly
  paramsInfo{end+1} = {'groupName',groupNames,'Name of group from which to do pRF analysis'};
  paramsInfo{end+1} = {'saveName','pRF','File name to try to save as'};
  paramsInfo{end+1} = {'restrict',restrict,'Restrict to the analysis to some subset of voxels. If you choose a base anatomy then it will restrict to the voxels that are on the base. If you choose an roi it will restrict the analysis to the voxels in the roi'};

  % if we give the option to continue an analysis
  if ~isempty(pRFAnalyses) && ~defaultParams
    continueParamsInfo{1} = {'continueAnalysis',0,'type=checkbox','Continue running a previously run analysis with same parameters. This is usually done on a different restriction set and will look at the old analysis to make sure not to compute voxels that have already been run. In the end will merge the analyses'};
    continueParamsInfo{2} = {'continueWhich',pRFAnalyses,'Which analysis to continue','contingent=continueAnalysis'};
    % add on restrict
    for i = 3:length(paramsInfo)
      continueParamsInfo{2+i-2} = paramsInfo{i};
    end
    for i = 3:length(continueParamsInfo)
      continueParamsInfo{i}{end+1} = 'contingent=continueAnalysis';
    end
    % put up dialog box with possibility to continue analysis
    continueParams = mrParamsDialog(continueParamsInfo,'Continue existing analysis?');
    if ~isempty(continueParams) && continueParams.continueAnalysis
      % get the parameters to continue from
      params = viewGet(v,'analysisParams',viewGet(v,'analysisNum',continueParams.continueWhich));
      % copy relevant ones (i.e. what new thing to restrict on)
      params.restrict =  continueParams.restrict;
      % tell pRF to set merge analysis instead of asking
      params.mergeAnalysis = true;
      % now get a list of all finished voxels
      a = viewGet(v,'analysis',viewGet(v,'analysisNum',continueParams.continueWhich));
      if isfield(a,'d')
	for i = 1:length(a.d)
	  if isfield(a.d{i},'linearCoords')
	    params.computedVoxels{i} = a.d{i}.linearCoords;
	  end
	end
      end
      return
    end
  end
end

%all of these parameters are for pRFFit
paramsInfo{end+1} = {'rfType',{'gaussian','gaussian-hdr'},'Type of pRF fit. Gaussian fits a gaussian with x,y,width as parameters to each voxel. gaussian-hdr fits also the hemodynamic response with the parameters of the hdr as below.'};
paramsInfo{end+1} = {'betaEachScan',false,'type=checkbox','Compute a separate beta weight (scaling) for each scan in the concanetation. This may be useful if there is some reason to believe that different scans have different magnitude responses, this will allow the fit to scale the magnitude for each scan'};
paramsInfo{end+1} = {'algorithm',{'nelder-mead','levenberg-marquardt'},'Which algorithm to use for optimization. Levenberg-marquardt seems to get stuck in local minimum, so the default is nelder-mead. However, levenberg-marquardt can set bounds for parameters, so may be better for when you are trying to fit the hdr along with the rf, since the hdr parameters can fly off to strange values.'};
paramsInfo{end+1} = {'defaultConstraints',1,'type=checkbox','Sets how to constrain the search (i.e. what are the allowed range of stimulus parameters). The default is to constrain so that the x,y of the RF has to be within the stimulus extents (other parameter constrains will print to the matlab window). If you click this off a dialog box will come up after the stimulus has been calculated from the stimfiles allowing you to specify the constraints on the parameters of the model. You may want to custom constrain the parameters if you know something about the RFs you are trying to model (like how big they are) to keep the nonlinear fits from finding unlikely parameter estimates. Note that nelder-mead is an unconstrained fit so this will not do anything.'};
paramsInfo{end+1} = {'prefitOnly',false,'type=checkbox','Check this if you want to ONLY do a prefit and not optimize further. The prefit computes a preset set of model parameters (x,y,rfHalfWidth) and picks the one that produces a mdoel with the highest correlation with the time series. You may want to do this to get a quick but accurate fit so that you can draw a set of ROIs for a full analysis'};
paramsInfo{end+1} = {'quickPrefit',false,'type=checkbox','Check this if you want to do a quick prefit - this samples fewer x,y and rfWidth points. It is faster (especially if coupled with prefitOnly for a fast check), but the optimization routines may be more likely to get trapped into local minima or have to search a long time for the minimum'};
paramsInfo{end+1} = {'verbose',true,'type=checkbox','Display verbose information during fits'};
paramsInfo{end+1} = {'yFlipStimulus',0,'type=checkbox','Flip the stimulus image in the y-dimension. Useful if the subject viewed a stimulus through a mirror which caused the stimulus to be upside down in the y-dimension'};
paramsInfo{end+1} = {'xFlipStimulus',0,'type=checkbox','Flip the stimulus image in the x-dimension. Useful if the subject viewed a stimulus which was flipped in the x-dimension'};
paramsInfo{end+1} = {'timeShiftStimulus',0,'incdec=[-1 1]','Time shift the stimulus, this is useful if the stimulus created is not correct and needs to be shifted in time (i.e. number of volumes)'};
if ~isempty(v)
  paramsInfo{end+1} = {'dispStim',0,'type=pushbutton','buttonString=Display stimulus','Display the stimulus for scan number: dispStimScan with the current parameters','callback',@pRFGUIDispStimulus,'passParams=1','callbackArg',v};
  paramsInfo{end+1} = {'dispStimScan',viewGet(v,'curScan'),'incdec=[-1 1]',sprintf('minmax=[1 %i]',viewGet(v,'nScans')),'round=1','Sets which scans stimulus will be displayed when you press Display stimulus button'};
end
paramsInfo{end+1} = {'timelag',1,'minmax=[0 inf]','incdec=[-0.5 0.5]','The timelag of the gamma function used to model the HDR. If using gaussian-hdr, this is just the initial value and the actual value will be fit.'};
paramsInfo{end+1} = {'tau',0.6,'minmax=[0 inf]','incdec=[-0.1 0.1]','The tau (width) of the gamma function used to model the HDR. If using gaussian-hdr, this is just the initial value and the actual value will be fit.'};
paramsInfo{end+1} = {'exponent',6,'minmax=[0 inf]','incdec=[-1 1]','The exponent of the gamma function used to model the HDR. This is always a fixed param.'};
paramsInfo{end+1} = {'diffOfGamma',false,'type=checkbox','Set to true if you want the HDR to be a difference of gamma functions - i.e. have a positive and a delayed negative component'};
paramsInfo{end+1} = {'amplitudeRatio',0.3,'minmax=[0 inf]','incdec=[-0.1 0.1]','Ratio of amplitude of 1st gamma to second gamma','contingent=diffOfGamma'};
paramsInfo{end+1} = {'timelag2',2,'minmax=[0 inf]','incdec=[-0.5 0.5]','Time lag of 2nd ggamma for when you are using a difference of gamma functions','contingent=diffOfGamma'};
paramsInfo{end+1} = {'tau2',1.2,'minmax=[0 inf]','incdec=[-0.1 0.1]','The tau (width) of the second gamma function.','contingent=diffOfGamma'};
paramsInfo{end+1} = {'exponent2',6,'minmax=[0 inf]','incdec=[-1 1]','The exponent of the 2nd gamma function.','contingent=diffOfGamma'};
paramsInfo{end+1} = {'dispHDR',0,'type=pushbutton','buttonString=Display HDR','Display the HDR with the current parameters','callback',@pRFGUIDispHDR,'passParams=1'};
paramsInfo{end+1} = {'saveStimImage',0,'type=checkbox','Save the stim image back to the stimfile. This is useful in that the next time the stim image will not have to be recomputed but can be directly read from the file (it will get saved as a variable called stimImage'};
paramsInfo{end+1} = {'recomputeStimImage',0,'type=checkbox','Even if there is an already computed stim image (see saveStimImage) above, this will force a recompute of the image. This is useful if there is an update to the code that creates the stim images and need to make sure that the stim image is recreated'};
paramsInfo{end+1} = {'applyFiltering',1,'type=checkbox','If set to 1 then applies the same filtering that concatenation does to the model. Does not do any filtering applied by averages. If this is not a concat then does nothing besides mean subtraction. If turned off, will still do mean substraction on model.'};
paramsInfo{end+1} = {'stimImageDiffTolerance',5,'minmax=[0 100]','incdec=[-1 1]','When averaging the stim images should be the same, but some times we are off by a frame here and there due to inconsequential timing inconsistenices. Set this to a small value, like 5 to ignore that percentage of frames of the stimulus that differ within an average. If this threshold is exceeded, the code will ask you if you want to continue - otherwise it will just print out to the buffer the number of frames that have the problem'};

% Get parameter values
if defaultParams
  params = mrParamsDefault(paramsInfo);
else
  params = mrParamsDialog(paramsInfo,'Set pRF parameters');
end

% if empty user hit cancel
if isempty(params)
  if deleteViewOnExit,deleteView(v);end
  return
end

% just getting pRFFItParams, so we are done
if pRFFitParamsOnly,return,end

% get scans
v = viewSet(v,'groupName',params.groupName);
if ~isempty(scanList)
  params.scanNum = scanList;
elseif defaultParams
  params.scanNum = 1:viewGet(v,'nScans');
else
  params.scanNum = selectScans(v);
end
if isempty(params.scanNum)
  params = [];
  if deleteViewOnExit,deleteView(v);end
  return
end

if deleteViewOnExit,deleteView(v);end

% if we go here, split out the params that get passed to pRFFit
pRFFitParams = false;
for i = 1:length(paramsInfo)
  % Everything after the rfType is a param for pRFFit
  if pRFFitParams || strcmp(paramsInfo{i}{1},'rfType')
    % move params into pRFFit field
    params.pRFFit.(paramsInfo{i}{1}) = params.(paramsInfo{i}{1});
    params = rmfield(params,paramsInfo{i}{1});
    % all the next fields will be moved as well
    pRFFitParams = true;
    % remove these from the paramsInfo field
    if strcmp(paramsInfo{i}{1},'rfType')
      params.paramInfo = {params.paramInfo{1:i-1}};
      params.pRFFit.paramInfo = {params.paramInfo{i:end}};
    end
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% just display parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = dispParams(params)

paramsInfo = {};
% grab the parameters that are indicated in paramsInfo
if isfield(params,'paramInfo')
  % get the paramsInfo
  topParamsInfo = params.paramInfo;
  % go through each one
  for i = 1:length(topParamsInfo)
    % if it exists in the params filed then add it
    if isfield(params,topParamsInfo{i}{1})
      % and it to paramInfo
      paramsInfo{end+1} = params.paramInfo{i};
      % add the value from params 
      paramsInfo{end}{2} = params.(topParamsInfo{i}{1});
      % make it non editable
      paramsInfo{end}{end+1} = 'editable=0';
    end
  end
end

% add the pRFFit
if isfield(params,'pRFFit')
  pRFFitFieldNames = fieldnames(params.pRFFit);
  for iField = 1:length(pRFFitFieldNames)
    fieldName = pRFFitFieldNames{iField};
    if ~any(strcmp(fieldName,{'paramInfo','dispHDR','dispStim'}))
      paramsInfo{end+1} = {fieldName,params.pRFFit.(fieldName),'editable=0'};
    end
  end
end
retval = mrParamsDialog(paramsInfo,'pRF parameters');
retval = [];

%%%%%%%%%%%%%%%%%%%%%%%
%    pRFGUIDispHDR    %
%%%%%%%%%%%%%%%%%%%%%%%
function retval = pRFGUIDispHDR(params)

retval = [];

% compute for 25 seconds
t = 0:0.1:25;

% get the first gamma function
hdr = thisGamma(t,1,params.timelag,0,params.tau,params.exponent);
titleStr = sprintf('(timelag: %s tau: %s exponent: %s)',mlrnum2str(params.timelag),mlrnum2str(params.tau),mlrnum2str(params.exponent));

% if difference of gamma subtract second gamma from this one
if params.diffOfGamma
  hdr = hdr - thisGamma(t,params.amplitudeRatio,params.timelag2,0,params.tau2,params.exponent2);
  titleStr = sprintf('%s - %s x (timelag2: %s tau2: %s exponent2: %s)',titleStr,mlrnum2str(params.amplitudeRatio),mlrnum2str(params.timelag2),mlrnum2str(params.tau2),mlrnum2str(params.exponent2));
end
hdr = hdr/max(hdr);

% display
mlrSmartfig('pRFGUIDispHDR','reuse');clf;
plot(t,hdr,'k.-');
title(titleStr);
xlabel('Time (sec)');
ylabel('Amplitude');


%%%%%%%%%%%%%%%%%%%
%%   thisGamma   %%
%%%%%%%%%%%%%%%%%%%
function gammafun = thisGamma(time,amplitude,timelag,offset,tau,exponent)

exponent = round(exponent);
% gamma function
gammafun = (((time-timelag)/tau).^(exponent-1).*exp(-(time-timelag)/tau))./(tau*factorial(exponent-1));

% negative values of time are set to zero,
% so that the function always starts at zero
gammafun(find((time-timelag) < 0)) = 0;

% normalize the amplitude
if (max(gammafun)-min(gammafun))~=0
  gammafun = (gammafun-min(gammafun)) ./ (max(gammafun)-min(gammafun));
end
gammafun = (amplitude*gammafun+offset);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    pRFGUIDispStimulus    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = pRFGUIDispStimulus(v,params)

retval = [];

stim = pRFFit(v,params.dispStimScan,[],[],[],'justGetStimImage=1','fitTypeParams',params);
if ~isempty(stim)
  % concatenate all stim images
  im = [];
  for i = 1:length(stim)
    im = cat(3,im,stim{i}.im);
  end
  % display using mlrVol
  disp(sprintf('(pRFGUI) Flipping image in y dimension so that appears in mlrVol the way it was presented'));
  im = mlrImageXform(im,'flipY');
  mlrVol(im,'imageOrientation=1');
end

