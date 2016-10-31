% mrInit.m
%
%      usage: mrInit(<sessionParams>,<groupParams>,<justGetParams=1>,<defaultParams=1>)
%         by: justin gardner
%       date: 06/09/08
%    purpose: Init the session variables. usually just call with no arguments
%             and it will put up dialogs to set variables:
%
%             mrInit
%
%             This is an alternative to using mrInitRet (it does the same thing, but with
%             a different GUI).
%
%             You can also just get parameters:
% 
%             [sessionParams groupParams] = mrInit([],[],'justGetParams=1')
%             or
%             [sessionParams groupParams] = mrInit([],[],'justGetParams=1','defaultParams=1')
%
%             and then call mrInit again to set the parameters (useful for scripting)
%             mrInit(sessionParams,groupParams);
%
%             You can use mrSetPref to set preferences for magnet/coil and pulseSequence names
%             that will come down as choices in the GUI
%
function [sessionParams groupParams] = mrInit(sessionParams,groupParams,varargin)

% check arguments
getArgs(varargin,{'justGetParams=0','defaultParams=0','makeReadme=1','magnet=[]','coil=[]','pulseSequence=[]','subject=[]','operator=[]','description=[]','stimfileMatchList=[]'});

minFramePeriod = .01;  %frame period in sec outside which the user is prompted
maxFramePeriod = 100;  % that something weird's goin on

% get session params
if ieNotDefined('sessionParams')
  % get some defaults
  if isfile('mrSession.mat')
    load mrSession
    magnet = putOnTopOfList(session.magnet,mrGetPref('magnet'));
    coil = putOnTopOfList(session.coil,mrGetPref('coil'));
    pulseSequence = putOnTopOfList(session.protocol,mrGetPref('pulseSequence'));
    description = session.description;
    subject = session.subject;
    operator = session.operator;
  else
    if isempty('magnet') 
      magnet = mrGetPref('magnet');
    else 
      magnet = putOnTopOfList(magnet,mrGetPref('magnet'));
    end
    if isempty('coil') 
      coil = mrGetPref('coil');
    else
      coil = putOnTopOfList(coil,mrGetPref('coil'));
    end
    if isempty('pulseSequence') 
      pulseSequence = mrGetPref('pulseSequence');
    else
      pulseSequence = putOnTopOfList(pulseSequence,mrGetPref('pulseSequence'));
    end
    if isempty('description') description = '';end
    if isempty('subject') subject = '';end
    if isempty('operator') operator = '';end
  end
  % setup params dialog
  paramsInfo = {};
  paramsInfo{end+1} = {'description',description,'type=string','Description of the session. Can be anything to help you remember what the session was'};
  paramsInfo{end+1} = {'subject',subject,'type=string','Subject ID. Use an identifier that does not break the subject confidentiality'};
  paramsInfo{end+1} = {'operator',operator,'type=string','Person who operated the scanner'};
  paramsInfo{end+1} = {'magnet',magnet,'type=popupmenu','Choose which magnet you scanned on'};
  paramsInfo{end+1} = {'coil',coil,'type=popupmenu','Choose which coil you used'};
  paramsInfo{end+1} = {'pulseSequence',pulseSequence,'Choose which pulse sequence you scanned with'};
  paramsInfo{end+1} = {'pulseSequenceText','','Optional: enter some text to describe or qualify the pulseSequence text. This will get appended to the pulseSequence name chosen above'};

  % get the parameters
  if defaultParams
    sessionParams = mrParamsDefault(paramsInfo);
  else
    sessionParams = mrParamsDialog(paramsInfo,'Initialize session for mrTools');
  end
  
  % check for user cancel
  if isempty(sessionParams)
    sessionParams = [];groupParams = [];
    return
  end
 
end

% get scan params
if ieNotDefined('groupParams')
  % check for already existing mrSession
  if isfile('mrSession.mat')
    load mrSession;
    nScans = length(groups(1).scanParams);
    for i = 1:nScans
      scanNames{i} = groups(1).scanParams(i).fileName;
      descriptions{i} = groups(1).scanParams(i).description;
      totalFrames{i} = groups(1).scanParams(i).totalFrames;
      nFrames{i} = groups(1).scanParams(i).nFrames;
      junkFrames{i} = groups(1).scanParams(i).junkFrames;
    end
  else
    % get info about scans that live in Raw/TSeries
    scanDirName = fullfile('Raw','TSeries');
    scanFilenames = mlrImageGetAllFilenames(scanDirName,'mustNotHaveDotInFilename=1');
    nScans=0;scanNames = {};descriptions = {};totalFrames = {};nFrames = {};junkFrames = {};
    for i = 1:length(scanFilenames);
      % read the nifti header
      imageHeader = mlrImageHeaderLoad(fullfile(scanDirName,scanFilenames{i}));
      if imageHeader.nDim == 4
	nScans = nScans+1;
	scanNames{end+1} = scanFilenames{i};
	descriptions{end+1} = '';
	totalFrames{end+1} = imageHeader.dim(4);
	nFrames{end+1} = imageHeader.dim(4);
	junkFrames{end+1} = 0;
      else
	disp(sprintf('(mrInit) File %s is not a 4D Nifti file',scanFilenames{i}));
      end
    end
  end
  
  % check to see if we got any good scans
  if nScans == 0
    disp(sprintf('(mrInit) Could not find any valid scans in Raw/TSeries'));
    sessionParams = [];groupParams = [];
    return
  end
  
  % setup params dialog
  paramsInfo = {};
  paramsInfo{end+1} = {'scanNum',1,'incdec=[-1 1]',sprintf('minmax=[1 %i]',nScans),'The scanNumber','editable=0'};
  paramsInfo{end+1} = {'name',scanNames,'group=scanNum','type=string','editable=0','Names of scans'};
  paramsInfo{end+1} = {'totalFrames',totalFrames,'group=scanNum','type=numeric','Number of frames in scan','editable=0'};
  paramsInfo{end+1} = {'description',descriptions,'group=scanNum','type=string','Description of scans'};
  paramsInfo{end+1} = {'nFrames',nFrames,'group=scanNum','type=numeric','Number of frames to use in scan'};
  paramsInfo{end+1} = {'junkFrames',junkFrames,'group=scanNum','type=numeric','Number of frames to junk at beginning of each scan','callback',@mrInitJunkFrames};
  
  % get the stim files
  [stimFileNames stimFileVols] = getStimFiles;

  if ~isempty(stimFileNames)
    if isempty(stimfileMatchList)
      % match the stimfiles with the scans
      stimFileMatch = matchStimFiles(stimFileNames,stimFileVols,totalFrames);
    else
      % passed in stimFileMatch, need to make a cell array of cell arrays
      for i = 1:length(totalFrames)
	% for each scan, check to see if we have a matching entry in stimfileMatchList
	matchNum = [];
	if length(stimfileMatchList) >= i
	  % then go look for the name in our stimFileNames array
	  matchNum = find(strncmp(stimfileMatchList{i},stimFileNames,length(stimfileMatchList{i})));
	end
	% if there was no match, just choose the top most stimfile (this will create duplicates, I suppose
	% but it should not really happen if a user passes in the list
	if isempty(matchNum), matchNum = 1;end
	% put the matched one on top of list 
	stimFileMatch{i} = putOnTopOfList(stimFileNames{matchNum},stimFileNames);
      end
    end
    paramsInfo{end+1} = {'stimFile',stimFileMatch,'group=scanNum','Stimfile to use for this scan'};
  end

  % button to copy parameters
  paramsInfo{end+1} = {'copy','','type=pushbutton','callback',@mrInitCopyParameters,'buttonString=Copy parameters','Copy parameters to selected scans','passParams=1'};
  
  % bring up dialog
  if defaultParams
    groupParams = mrParamsDefault(paramsInfo);
  else
    groupParams = mrParamsDialog(paramsInfo,'Choose group parameters');
  end
  
  % check for user cancel
  if isempty(groupParams)
    sessionParams = [];groupParams = [];
    return
  end
end

if ~justGetParams
  if ~ieNotDefined('sessionParams') && ~ieNotDefined('groupParams')
    disp(sprintf('(mrInit) Saving session information'));
    weirdFramePeriods = zeros(1,length(groupParams.totalFrames));
    % create session variable
    session.mrLoadRetVersion = mrLoadRetVersion;
    session.description = sessionParams.description;
    session.subject = sessionParams.subject;
    session.operator = sessionParams.operator;
    session.magnet = sessionParams.magnet;
    session.coil = sessionParams.coil;
    session.protocol = sprintf('%s: %s',sessionParams.pulseSequence,sessionParams.pulseSequenceText);
    % create groups variables
    groups(1).name = 'Raw';
    scanParams = [];
    tseriesDir = 'Raw/TSeries';
    for iScan=1:length(groupParams.totalFrames)
      name = fullfile(tseriesDir, groupParams.name{iScan});
      hdr = mlrImageReadNiftiHeader(name);
      scanParams(iScan).dataSize = hdr.dim([2,3,4])';
      scanParams(iScan).description = groupParams.description{iScan};
      scanParams(iScan).fileName = groupParams.name{iScan};
      scanParams(iScan).originalFileName{1} = groupParams.name{iScan}; % otherwise motionComp has problems 
      scanParams(iScan).originalGroupName{1} = groups(1).name; % ditto
      scanParams(iScan).fileType = 'Nifti';
      niftiSpaceUnit = bitand(hdr.xyzt_units, hex2dec('07')); 
      niftiTimeUnit = bitand(hdr.xyzt_units,hex2dec('38'));
      if niftiTimeUnit == 8 % seconds
	scanParams(iScan).framePeriod = hdr.pixdim(5)./1;
      elseif niftiTimeUnit == 16 % milliseconds
	scanParams(iScan).framePeriod = hdr.pixdim(5)./1000;
      elseif niftiTimeUnit == 24 % microseconds
	scanParams(iScan).framePeriod = hdr.pixdim(5)./10e6;
      end
      %detect weird frame periods
      if scanParams(iScan).framePeriod>maxFramePeriod 
        weirdFramePeriods(iScan) = scanParams(iScan).framePeriod;
      elseif scanParams(iScan).framePeriod<minFramePeriod
        weirdFramePeriods(iScan) = scanParams(iScan).framePeriod;
      end
      scanParams(iScan).junkFrames = groupParams.junkFrames(iScan);
      scanParams(iScan).nFrames = groupParams.nFrames(iScan);
      scanParams(iScan).niftiHdr = hdr;
      scanParams(iScan).totalFrames = hdr.dim(5);
      scanParams(iScan).voxelSize = hdr.pixdim([2,3,4])';
    end  
    
    %if weird frame periods have been detected, ask user
    if any(weirdFramePeriods)
      scanParams = fixFramePeriods(scanParams,weirdFramePeriods,minFramePeriod,maxFramePeriod,tseriesDir);
    end
    
    % put into group variable
    for iScan=1:length(groupParams.totalFrames)
      [tf groups(1).scanParams(iScan)] = isscan(scanParams(iScan));
      if ~tf
	disp(sprintf('(mrInit) Scan params for %i did not validate',iScan));
      end
      % check for stimfiles field and add to aux params if found
      if isfield(groupParams,'stimFile')
	if ~strcmp(groupParams.stimFile{iScan},'None')
	  % remove descriptive text when saving
	  groups(1).auxParams(iScan).stimFileName = strtok(groupParams.stimFile{iScan},':');
	else
	  groups(1).auxParams(iScan).stimFileName = '';
	end	  
      end
    end
    
    % tag on auxParams to groups if it wasn't set above
    if ~isfield(groups(1),'auxParams')
      groups(1).auxParams = [];
    end
    
    % check for mrSession
    if isfile('mrSession.mat')
      if askuser('(mrInit) mrSession.mat already exists. Overwrite?');
	disp(sprintf('(mrInit) Copying old mrSession.mat mrSession.old.mat'));
	movefile('mrSession.mat','mrSession.old.mat');
	disp(sprintf('(mrInit) Saving new mrSession'));
	save mrSession session groups;
	% disp(sprintf('(mrInit) Creating new Readme'));
	if makeReadme
	  mrReadme(session, groups);
	end
      end
    else
      disp(sprintf('(mrInit) Saving new mrSession'));
      save mrSession session groups;
      disp(sprintf('(mrInit) Creating Readme'));
      if makeReadme
	mrReadme(session, groups);
      end
    end
  else
    disp(sprintf('(mrInit) No session information saved'));
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   mrInitCopyParameters  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function params = mrInitCopyParameters(params)

% put a button dialog to choose which scans to copy parameters to
copyTo = buttondlg('Copy parameters to scans:',params.name);

% user cancel
if isempty(copyTo),return,end

% go copy parameters
for i = 1:length(copyTo)
  if copyTo(i)
    params.description{i} = params.description{params.scanNum};
    params.junkFrames(i) = params.junkFrames(params.scanNum);
    params.nFrames(i) = params.totalFrames(i)-params.junkFrames(i);
  end
end
mrParamsSet(params);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mrInitJunkFrames    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function params = mrInitJunkFrames(params)

% reset the nFrames field
params.nFrames(params.scanNum) = params.totalFrames(params.scanNum)-params.junkFrames(params.scanNum);
mrParamsSet(params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% match stimfiles with scans  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimFileMatch = matchStimFiles(stimFileNames,stimFileVols,totalFrames)

stimFileMatch = {};
% for now, just match in order
for i = 1:length(stimFileNames)
  % make sure there are enough vols (just check to see if there
  % are greater than 1/10 the total frames) since with tSense
  % etc there may be more totalFrames than vols, so this is
  % just to reject stimfiles with only a few volumes
  if (length(stimFileVols) < i) || (length(totalFrames)<i) || (stimFileVols(i) > totalFrames{i}/10)
    stimFileMatch{end+1} = putOnTopOfList(stimFileNames{i},{stimFileNames{:},'None'});
    % only make enough matches for each scan
    if length(stimFileMatch)==length(totalFrames)
      break;
    end
  end
end

% all unmatched scans
for i = length(stimFileMatch)+1:length(totalFrames)
  stimFileMatch{end+1} = putOnTopOfList('None',cellcat(stimFileNames,{'None'}));
end

%%%%%%%%%%%%%%%%%%%%%%
%  get stimFileNums  %
%%%%%%%%%%%%%%%%%%%%%%
function [stimFileNames stimFileVols] = getStimFiles

stimfileDir = 'Etc';
dirList = dir(stimfileDir);
stimFileNums = [];stimFileNames = {};stimFileVols = [];
for i = 1:length(dirList);
  name = dirList(i).name;
  % name should be of form yymmdd_stimnn.mat (nn is sequence number)
  if (regexp(name,'\d\d\d\d\d\d_\w\w\w\w\d\d.mat'))
    % get the stimfile number
    stimFileNums(end+1) = str2num(name(12:13));
    % try to load and get some info
    load(sprintf(fullfile(stimfileDir,name)));
    if ieNotDefined('myscreen')
      disp(sprintf('(mrInit) %s does not look like an mgl stimfile (does not have a myscreen variable)',name));
      continue;
    end
    % make a string of some info myscreen
    stimfileStr = '';
    if isfield(myscreen,'volnum')
      stimfileStr = sprintf('%s[%i vols] ',stimfileStr,myscreen.volnum);
      stimFileVols(end+1) = myscreen.volnum;
    else
      stimFileVols(end+1) = nan;
    end
    if isfield(myscreen,'starttime')
      stimfileStr = sprintf('%s%s ',stimfileStr,myscreen.starttime);
    end
    if isfield(myscreen,'endtime')
      stimfileStr = sprintf('%s(End: %s) ',stimfileStr,myscreen.endtime);
    end
    stimFileNames{end+1} = sprintf('%s: %s',name,stimfileStr);
  end
end
