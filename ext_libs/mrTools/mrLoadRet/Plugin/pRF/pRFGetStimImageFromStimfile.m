% pRFGetStimImageFromStimfile.m
%
%        $Id:$ 
%      usage: stim = pRFGetStimImageFromStimfile(stimfile,<timePoints>)
%         by: justin gardner
%       date: 10/17/11
%    purpose: Pass in a stimfile (can be either a string filename, or a strucutre
%             with myscreen/task) created by mglRetinotopy (make sure this is a
%             stimfile created with a version of mglRetinotopy past 10/2011 which
%             has the proper variables stored to enable reconstruction). Will
%             create a volume of dimensions x,y,t with the stimulus image (load
%             in mlrVol to view). stim.x and stim.y are the X and Y coordinates
%             in degrees of every point. stim.t is the array of times at
%             which image is taken.
%
%             Optionally arguments:
%
%             timePoints: array for which the stim image should be computed. 
%
%             Note for developers - this function needs to keep up-to-date with
%             any changes in the display loop of mglRetinotopy to interpret
%             the stimfiles correctly
%
function stim = pRFGetStimImageFromStimfile(stimfile,varargin)

% set default return arguments
stim = [];

% check arguments
if nargin < 1
  help pRFGetStimImageFromStimfile
  return
end

% parse arguments
timePoints = [];screenWidth = [];screenHeight = [];volTrigRatio = [];
xFlip = [];yFlip = [];timeShift = [];verbose = [];
getArgs(varargin,{'timePoints=[]','screenWidth=[]','screenHeight=[]','volTrigRatio=[]','xFlip=0','yFlip=0','timeShift=0','verbose=1','saveStimImage=0','recomputeStimImage=0'});

% handle cell array
if iscell(stimfile) && ((length(stimfile)>1) || (length(stimfile{1})>1))
  for i = 1:length(stimfile)
    % get current volTrigRatio
    if isempty(volTrigRatio)
      thisVolTrigRatio = [];
    else
      thisVolTrigRatio = volTrigRatio{i};
    end
    stim{i} = pRFGetStimImageFromStimfile(stimfile{i},'timePoints',timePoints,'screenWidth',screenWidth,'screenHeight',screenHeight,'volTrigRatio',thisVolTrigRatio,'xFlip',xFlip,'yFlip',yFlip,'timeShift',timeShift,'verbose',verbose,'saveStimImage',saveStimImage,'recomputeStimImage',recomputeStimImage);
    if isempty(stim{i}),stim = [];return;end
  end
  return
end

% check volTrigRatio
if iscell(volTrigRatio)
  if length(volTrigRatio) > 1
    disp(sprintf('(pRFGetStimImageFromStimfile) volTrigRatio should not be of length greater than one (length=%i) using only the first value of %i',length(volTrigRatio),volTrigRatio{1}));
  end
  volTrigRatio = volTrigRatio{1};
end

% get the stimfile
s = getStimfile(stimfile);
if isempty(s),return,end

% check that we have a stimfile that is interpretable
% by this program
[tf s taskNum] = checkStimfile(s);
if ~tf,return,end

% check to see if a stimImage exists
if ~isfield(s,'pRFStimImage') || recomputeStimImage
  % make the traces
  s.myscreen = makeTraces(s.myscreen,verbose);

  % get task variables
  e = getTaskParameters(s.myscreen,s.task{taskNum});

  % get some traces of things of interest
  s.time = s.myscreen.time;
  s.maskPhase = s.myscreen.traces(s.task{taskNum}{1}.maskPhaseTrace,:);
  s.blank = s.myscreen.traces(s.task{taskNum}{1}.blankTrace,:);
  s.vol = s.myscreen.traces(1,:);
  s.trialVol = e.trialVolume;
  s.trialTicknum = e.trialTicknum;
  s.blank = e.randVars.blank;
  if any(s.stimulus.stimulusType == [3 4])
    s.barAngle = e.parameter.barAngle;
    s.elementAngle = e.randVars.elementAngle;
  end

  % if no timepoints, then get one for each volume
  if isempty(timePoints)
    timePoints = s.time(find(s.vol))-s.time(first(find(s.vol)));
    
    % add timepoints for sense acceleration
    if ~isempty(volTrigRatio)
      framePeriod = median(diff(timePoints));
      accTimePoints = [];
      for i = 1:length(timePoints);
	accTimePoints(end+1) = timePoints(i);
	for j = 1:(volTrigRatio-1)
	  accTimePoints(end+1) = timePoints(i)+framePeriod*j/volTrigRatio;
	end
      end
      timePoints = accTimePoints;
    end
    % add one frame extra since the stimulus usually changes after
    % the volume
    timePoints = timePoints + s.myscreen.frametime;
  end
  stim.t = timePoints;

  % check the pixel dimensions
  if ~isempty(screenWidth) && ~isempty(screenHeight)
  elseif isempty(screenWidth) && isempty(screenHeight)
    screenWidth = round(s.myscreen.screenWidth/10);
    screenHeight = round(s.myscreen.screenHeight/10);
  elseif isempty(screenWidth)
    % no screenWith, compute based on aspect ratio
    screenWidth = round(screenHeight*s.myscreen.screenWidth/s.myscreen.screenHeight);
  elseif isempty(screenHeight)
    screenHeight = round(screenWidth*s.myscreen.screenHeight/s.myscreen.screenWidth);
  end  

  % open the screen
  if mglGetParam('displayNumber') ~= -1,mglClose;end
  mglSetParam('offscreenContext',1);
  mglOpen(0,screenWidth,screenHeight);
  mglVisualAngleCoordinates(s.myscreen.displayDistance,s.myscreen.displaySize);

  % create stim.x and stim.y
  imageWidth = s.myscreen.imageWidth;
  imageHeight = s.myscreen.imageHeight;
  [stim.x stim.y] = ndgrid(-imageWidth/2:imageWidth/(screenWidth-1):imageWidth/2,-imageHeight/2:imageHeight/(screenHeight-1):imageHeight/2);

  if verbose,disppercent(-inf,'(pRFGetStimImageFromStimfile) Computing stimulus images');end
  warnOnStimfileMissingInfo = true;
  for iImage = 1:length(stim.t)
    im = createMaskImage(s,stim.t(iImage));
    % if no image, that probably means the stimfile ended early
    if isempty(im)
      % check to see if we are within (arbitrarily) 5% of the end
      % and use that as a cutoff for asking the user if something drastically
      % wrong has occurred.
      if warnOnStimfileMissingInfo
	if (1-iImage/length(stim.t)) > 0.05
	  if askuser('Your stimfile is missing information for volume %i of %i. This might be because you have linked the wrong stimfile or that the stim program ended before the scan or that there is some other problem with the stimfile. It would be a good idea to try to fix this cause this may now be generating the wrong stimulus. Continue anyway?',0,1)
	    warnOnStimfileMissingInfo = false;
	  else
	    % user did not agree to continue, bail out
            stim=[];
	    return
	  end
	end
      end
      % just put up the warning
      disp(sprintf('(pRFGetStimImageFromStimfile) !!! Missing stimulus info for volume %i of %i. Setting to blank image. !!!',iImage,length(stim.t)));
      im = zeros(screenWidth,screenHeight);
    end
    stim.im(1:screenWidth,1:screenHeight,iImage) = im;
    if verbose,disppercent(iImage/length(stim.t));end
  end
  if verbose,disppercent(inf);end

  % close screen
  mglSetParam('offscreenContext',0);
  mglClose;
  % save stim back to stimFile if called for
  if saveStimImage
    saveStimImageToStimfile(stim,stimfile);
  end
else
  % stim image was stored, just reclaim it
  disp(sprintf('(pRFGetStimImageFromStimfile) Loaded stim image from stimfile.'));
  stim = s.pRFStimImage;
end

% apply transformations if called for
if xFlip
  disp(sprintf('(pRFGetStimImageFromStimfile) X flipping stimulus image'));
  for i = 1:size(stim.im,3)
    stim.im(:,:,i) = flipud(stim.im(:,:,i));
  end
end
if yFlip
  disp(sprintf('(pRFGetStimImageFromStimfile) Y flipping stimulus image'));
else
  % note that we actually flip here by default since that seems to make
  % everything right in the analysis.
  for i = 1:size(stim.im,3)
    stim.im(:,:,i) = fliplr(stim.im(:,:,i));
  end
end
if timeShift
  disp(sprintf('(pRFGetStimImageFromStimfile) Time shifting stimulus image by %i',timeShift));
  stim.im = circshift(stim.im,[0 0 timeShift]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%    createMaskImage    %
%%%%%%%%%%%%%%%%%%%%%%%%%
function maskImage = createMaskImage(s,t)

% find the beginning of the experiment
firstTimepoint = find(s.vol);
firstTimepoint = firstTimepoint(1);

% find the timepoint that corresponds to this time
thisTimepoint = s.time(firstTimepoint)+t;
thisTimepoint = find(thisTimepoint <= s.time);
if isempty(thisTimepoint)
  disp(sprintf('(pRFGetStimImageFromStimfile) Timepoint %0.1fs does not exist in stimfile. This might have happened if you have linked the wrong stimfile with the scan - in which case, setStimfile to change stimfile linking',t));
  maskImage = [];
  return
end
% get timepoint
thisTimepoint = thisTimepoint(1);

% get current volume number
thisVol = cumsum(s.vol);
thisVol = thisVol(thisTimepoint);

% get curent trial number
thisTrial = find(s.trialVol <= thisVol);
if isempty(thisTrial)
  disp(sprintf('(pRFGetStimImageFromStimfile) Trial for volume %i not found. The first trial starts at volume %i. Consider discarding the scan associated with stimfile %s',thisVol,first(s.trialVol),getLastDir(s.myscreen.stimfile)));
  maskImage = [];
  return;
end
thisTrial = thisTrial(end);

% make sure that the timepoint is valid for the trial
thisTimepoint = max(thisTimepoint,s.trialTicknum(thisTrial));
if length(s.trialTicknum) > thisTrial
  thisTimepoint = min(thisTimepoint,s.trialTicknum(thisTrial+1));
end

% pull out stimulus variable
global stimulus;
stimulus = s.stimulus;

if isfield(s,'barAngle')
  % now make a rotation matrix for the background angle
  elementAngle = s.elementAngle(thisTrial);
  co = cos(pi*elementAngle/180);
  si = sin(pi*elementAngle/180);
  stimulus.elementRotMatrix = [co si;-si co];

  % now make a rotation matrix for the bar angle we want to present
  barAngle = s.barAngle(thisTrial);
  co = cos(pi*barAngle/180);
  si = sin(pi*barAngle/180);
  stimulus.maskBarRotMatrix = [co si;-si co];
  
  % see whether this is a blank
  if barAngle == -1,blank = true;else blank = false;end

  % clear screen
  mglClearScreen(1);
else
  % for rings and wedges see if it is a blank
  blank = s.blank(thisTrial);
  if blank
    disp(sprintf('(pRFGetStimImageFromStimfile) Blank trials not yet coded here'));
    % this just needs to read the segment and decide which half
    % of the trial to balnk out
    keyboard
  end
  mglClearScreen(0.5);
  mglFillOval(0,0,[stimulus.maxRadius*2 stimulus.maxRadius*2],[1 1 1]);
end

% set the current mask
stimulus.currentMask = s.maskPhase(thisTimepoint);
if stimulus.currentMask == 0,blank = true;end

% update bars only if this is not a blank frame
if ~blank
  updateRetinotopyStimulus(stimulus,s.myscreen);
else
  % clear screen
  mglClearScreen(0.5);
end

% grab the screen
maskImage = mglFrameGrab;

% make into a black and white image
maskImage((maskImage > 0.51) | (maskImage < 0.49)) = 1;
maskImage((maskImage < 0.51) & (maskImage > 0.49)) = 0;
maskImage = maskImage(:,:,1);

% DEBUG CODE - will draw each frame of the stimulus to a figure

%disp(sprintf('Trial %i maskPhase: %i',thisTrial,stimulus.currentMask))
%mlrSmartfig('pRF','reuse');clf;imagesc(maskImage);
%title(sprintf('Trial %i maskPhase: %i',thisTrial,stimulus.currentMask))
%drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to draw retinotopy stimulus to screen
% this function has been taken out of mglRetinotopy
% the only thing changed is that the mglQuad call has been
% commented out so that the background pattern 
% does not draw
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = updateRetinotopyStimulus(stimulus,myscreen)

if any(stimulus.stimulusType == [3 4])

  % update the phase of the sliding wedges
  stimulus.phaseNumRect = 1+mod(stimulus.phaseNumRect,stimulus.nRect);

  % draw the whole stimulus pattern, rotate to the element angle
  x = stimulus.xRect{stimulus.phaseNumRect};
  y = stimulus.yRect{stimulus.phaseNumRect};
  coords(1:2,:) = stimulus.elementRotMatrix*[x(1,:);y(1,:)];
  coords(3:4,:) = stimulus.elementRotMatrix*[x(2,:);y(2,:)];
  coords(5:6,:) = stimulus.elementRotMatrix*[x(3,:);y(3,:)];
  coords(7:8,:) = stimulus.elementRotMatrix*[x(4,:);y(4,:)];
%  mglQuad(coords(1:2:8,:),coords(2:2:8,:),stimulus.cRect{stimulus.phaseNumRect},1);

  % compute the center of the bar
  barCenter = repmat(stimulus.barCenter(stimulus.currentMask,:),size(stimulus.maskBarLeft,1),1);
  % compute the left and right masks (covering up everything except the bar)
  % by shifting by the barCenter and rotating the coordinates for the angle we want
  maskBarLeft = stimulus.maskBarRotMatrix*(barCenter+stimulus.maskBarLeft)';
  maskBarRight = stimulus.maskBarRotMatrix*(barCenter+stimulus.maskBarRight)';

  % draw the bar masks
  mglPolygon(maskBarLeft(1,:),maskBarLeft(2,:),0.5);
  mglPolygon(maskBarRight(1,:),maskBarRight(2,:),0.5);
else
  % update the phase of the sliding wedges
  stimulus.phaseNum = 1+mod(stimulus.phaseNum,stimulus.n);
  % draw the whole stimulus pattern
%  mglQuad(stimulus.x{stimulus.phaseNum},stimulus.y{stimulus.phaseNum},stimulus.c{stimulus.phaseNum},1);
  
  % mask out to get a wedge
  if stimulus.stimulusType == 1
    mglPolygon(stimulus.maskWedgeX{stimulus.currentMask},stimulus.maskWedgeY{stimulus.currentMask},0.5);
    % or mask out to get a ring
  else
    mglPolygon(stimulus.maskInnerX{stimulus.currentMask},stimulus.maskInnerY{stimulus.currentMask},0.5);
    mglQuad(stimulus.maskOuterX{stimulus.currentMask},stimulus.maskOuterY{stimulus.currentMask},stimulus.maskOuterC{stimulus.currentMask});
  end
end

%%%%%%%%%%%%%%%%%%%%%
%    getStimfile    %
%%%%%%%%%%%%%%%%%%%%%
function s = getStimfile(stimfile)

s = [];

% deal with a cell array of stimfiles (like in an average)
if iscell(stimfile)
  for i = 1:length(stimfile)
    s{i} = getStimfile(stimfile{i});
    if isempty(s{i}),return;end
  end
  return
end

% load stimfile
if isstr(stimfile)
  stimfile = setext(stimfile,'mat');
  if ~isfile(stimfile)
    disp(sprintf('(pRFGetStimImageFromStimfile) Could not open stimfile: %s',stimfile));
    return
  end
  s = load(stimfile);
elseif isstruct(stimfile)
  % see if this is a myscreen
  if isfield(stimfile,'imageWidth')
    % check for task field
    if isfield(stimfile,'task')
      s.task = stimfile.task;
      stimfile = rmfield(stimfile,'task');
    end
    % check for stimulus field
    if isfield(stimfile,'stimulus')
      s.stimulus = stimfile.stimulus;
      stimfile = rmfield(stimfile,'stimulus');
    end
    % set myscreen field
    s.myscreen = stimfile;
  % else a variable with myscreen, task and stimulus or pRFStimImage
  elseif isfield(stimfile,'myscreen') || isfield(stimfile,'pRFStimImage')
    % copy fields over
    if isfield(stimfile,'myscreen')
      s.myscreen = stimfile.myscreen;
    end
    if isfield(stimfile,'task')
      s.task = stimfile.task;
    end
    if isfield(stimfile,'stimulus')
      s.stimulus = stimfile.stimulus;
    end
    if isfield(stimfile,'pRFStimImage')
      s.pRFStimImage = stimfile.pRFStimImage;
    end
  end
end

% if you have a pRFStimImage then don't bother with the rest of the fields
if ~isfield(s,'pRFStimImage')
  % check fields
  checkFields = {'myscreen','task','stimulus'};
  for i = 1:length(checkFields)
    if ~isfield(s,checkFields{i})
      stimfileName = '';
      if isfield(s,'myscreen') && isfield(s.myscreen,'stimfile')
	stimfileName = getLastDir(s.myscreen.stimfile);
      end
      disp(sprintf('(pRFGetStimImageFromStimfile) !!! Missing variable: %s in stimfile %s !!!',checkFields{i},stimfileName));
      s = [];
      return
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%
%    checkStimfile    %
%%%%%%%%%%%%%%%%%%%%%%%
function [tf s taskNum] = checkStimfile(s)

tf = true;
s = cellArray(s);
taskNum = [];

stimulusType = [];
barAngle = [];
direction = [];

for i = 1:length(s)
  thiss = s{i};
  if isempty(thiss)
    disp(sprintf('(pRFGetStimImageFromStimfile) Missing stimfile'));
    tf = false;
    return
  end
  % if this has a pRFStimImage then we are ok
  if isfield(thiss,'pRFStimImage')
    continue;
  end
  dispstr = sprintf('%s: vols=%i',thiss.myscreen.stimfile,thiss.myscreen.volnum);
  % first check if this is a retinotpy stimfile - it should
  % have a task which is mglRetinotopy
  taskNum = [];
  for iTask = 1:2
    if (length(thiss.task) >= iTask) && (isequal(thiss.task{iTask}{1}.taskFilename,'mglRetinotopy.m') || isequal(thiss.task{iTask}{1}.taskFilename,'gruRetinotopy.m'))
      taskNum = iTask;
    end
  end
  if isempty(taskNum)
    disp(sprintf('(pRFGetStimImageFromStimfile:checkStimfile) Stimfile: %s',dispstr));
    disp(sprintf('(pRFGetStimImageFromStimfile:checkStimfile) The stimfile does not appear to have been created by mglRetinotopy'));
    tf = false;
    return
  end

  % check for proper saved fields
  missing = '';
  if ~isfield(thiss.task{taskNum}{1},'randVars') missing = 'randVars';end
  if ~isfield(thiss.task{taskNum}{1},'parameter') missing = 'parameter';end
  if ~any(strcmp('maskPhase',thiss.myscreen.traceNames)) missing = 'maskPhase';end
  if ~any(strcmp('blank',thiss.myscreen.traceNames)) missing = 'blank';end
  if ~isempty(missing)
    disp(sprintf('(pRFGetStimImageFromStimfile:checkStimfile) Stimfile: %s',dispstr));
    disp(sprintf('(pRFGetStimImageFromStimfile:checkStimfile) The stimfile does not appear to have been created by the latest version of mglRetinotopy which contains the field %s necessary for reconstructing the stimulus. Consider running a dummy run with a newer version of mglRetinotpy with the same parameters (see mglSimulateRun to simulate backticks) and then use that stimfile instead of this one.',missing));
    tf = false;
    return
  end

  % check for necessary variables
  e = getTaskParameters(thiss.myscreen,thiss.task{taskNum}{1});

  % now check for each variable that we need
  varnames = {'blank'};
  for i = 1:length(varnames)
    varval = getVarFromParameters(varnames{i},e);
    if isempty(varval)
      disp(sprintf('(pRFGetStimImageFromStimfile:checkStimfile) Stimfile: %s',dispstr));
      disp(sprintf('(pRFGetStimImageFromStimfile:checkStimfile) The stimfile does not appear to have been created by the latest version of mglRetinotopy which contains the variable %s necessary for reconstructing the stimulus. Consider running a dummy run with a newer version of mglRetinotpy with the same parameters (see mglSimulateRun to simulate backticks) and then use that stimfile instead of this one',varnames{i}));
      tf = false;
      return
    end
  end
  
  % check for matching stimfiles
  if ~isempty(stimulusType) && (stimulusType ~= thiss.stimulusType)
    disp(sprintf('(pRFGetStimImageFromStimfile:checkStimfile) !!! Stimfile %s does not match previous one !!! Have you averaged together scans with different stimulus conditions?'));
  end
  if any(thiss.stimulus.stimulusType == [3 4])
    varval = getVarFromParameters('barAngle',e);
    if ~isempty(barAngle) && ~isequal(varval,barAngle)
      disp(sprintf('(pRFGetStimImageFromStimfile:checkStimfile) !!! Stimfile %s does not match previous one !!! The barAngles are different! Have you averaged together scans with different stimulus conditions?'));
    end
    barAngle = varval;
  else
    if ~isempty(direction) && (thiss.stimulus.direction ~= direction)
      disp(sprintf('(pRFGetStimImageFromStimfile:checkStimfile) !!! Stimfile %s does not match previous one !!! The directions are different! Have you averaged together scans with different stimulus conditions?'));
    end
    direction = thiss.stimulus.direction;
  end
end

s = s{end};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    saveStimImageToStimfile    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveStimImageToStimfile(stim,stimfile)

% make sure stimfile is a cell array
stimfile = cellArray(stimfile);

% first reload the stimfile
for iStimfile = 1:length(stimfile)
  if isfield(stimfile{iStimfile},'filename')
    s = load(stimfile{iStimfile}.filename);
    if isempty(s)
      disp(sprintf('(pRFGetStimImageFromStimfile:saveStimImageToStimfile) Could not load stimfile %s. Unable to save stim image back to stimfile',stimfile{iStimfile}.filename));
    else
      % append the stim image and save back
      s.pRFStimImage = stim;
      save(stimfile{iStimfile}.filename,'-struct','s');
      disp(sprintf('(pRFGetStimImageFromStimfile:saveStimImageToStimfile) Saved stimImage to %s.',stimfile{iStimfile}.filename));
    end
  else
    disp(sprintf('(pRFGetStimImageFromStimfile:saveStimImageToStimfile) Missing filename in stimfile structure, could not save stimImage back to stimfile'));
  end
end

