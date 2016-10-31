function [val val2] = viewGet(view,param,varargin)
%
%   val = viewGet(view,param,varargin)
%
% Read the parameters of a view
% Access to these structures should go through this routine and through
% viewSet.
%
% ====================================================================
% Type viewGet w/out any arguments for the full set of parameters and for
% optional arguments. If you want help on a specific parameter, you
% can also do:
% viewGet([],'parameterName','help')
% ====================================================================
%
% view can be either a view structure or a viewNum which is interpreted as
% MLR.views{viewNum}. For some params (e.g., 'subject', 'groupNum'), view
% can be [].
%
% For some params, you must pass additional arguments in varargin. For
% some params, the additional arguments are optional. An example that
% illustrates all of these possibilities is:
%       n = viewGet(view,'nScans',[groupNum]);
% where view can be either view structure or []. If view is a view
% structure, then and groupNum is optional (defaults to the current group).
%
% Examples:
%
% tseriesdir = viewGet(view,'tseriesdir');
% n = viewGet(view,'numberofGroups');
% groupNames = viewGet([],'groupNames');
% groupNum = viewGet(view,'currentGroup');
% n = viewGet(view,'nScans',groupNum);
% n = viewGet([],'nScans',groupNum);
%
% 6/2004 djh
% 11/2006 jlg added help and
% originalScanNum,originalFileName,originalGroupName, scanNum,
% dicomName, dicom, stimFileName, stimFile, concatInfo, transforms, TR
%	$Id: viewGet.m 2868 2013-09-27 04:17:24Z justin $

mrGlobals

% if ieNotDefined('view'), error('No view defined.'); end
if ieNotDefined('param')
  dispViewGetHelp;
  return
elseif ~ieNotDefined('varargin') && (length(varargin)==1) && isstr(varargin{1}) && (strcmp(lower(varargin{1}),'help') || strcmp(lower(varargin{1}),'?'))
  dispViewGetHelp(param);
  return
end
% If 'view' is a viewNum, find the corresponding view structure
if isnumeric(view) && isscalar(view) && ~isempty(view) && (view<= length(MLR.views))
  view = MLR.views{view};
end

% Initialize return value
val = [];val2 = [];

switch lower(param)
  
  case {'view'}
    % view = viewGet(view,'view')
    % view = viewGet([],'view',viewNum)
    if length(varargin) == 1
      if (varargin{1} > 0) & (varargin{1} <= length(MLR.views))
        val = MLR.views{varargin{1}};
      end
    else
      val = view;
    end
  case {'viewtype','type'}
    % viewType = viewGet(view,'viewType')
    val = view.viewType;
  case {'viewnums'}
    % Returns an array containing the viewNum of every active view.
    % Returns empty if there are no open views
    % viewNums = viewGet([],'viewnums')
    for v=1:length(MLR.views)
      if isview(MLR.views{v})
        val = [val,v];
      end
    end
  case {'viewnum'}
    % viewNum = viewGet(view,'viewNum')
    val = view.viewNum;
  case {'subject'}
    % subject = viewGet(view,'subject')
    val = MLR.session.subject;
  case {'sessiondescription'}
    % subject = viewGet(view,'sessionDescription')
    val = MLR.session.description;
  case {'operator'}
    % subject = viewGet(view,'operator')
    val = MLR.session.operator;
  case {'magnet'}
    % subject = viewGet(view,'magnet')
    val = MLR.session.magnet;
  case {'coil'}
    % subject = viewGet(view,'coil')
    val = MLR.session.coil;
  case {'protocol'}
    % subject = viewGet(view,'protocol')
    val = MLR.session.protocol;
    
    % subdirectories
  case {'homedir','homedirectory','sessiondirectory'}
    % homeDir = viewGet(view,'homeDir')
    val = MLR.homeDir;
  case {'anatomydir'}
    % anatomyDir = viewGet(view,'anatomydir')
    homeDir = viewGet(view,'homedir');
    val = fullfile(homeDir,'Anatomy');
    if ~exist(val,'dir')
      mkdir(homeDir,'Anatomy');
    end
  case {'datadir'}
    % dataDir = viewGet(view,'datadir',[groupNum])
    % dataDir = viewGet(view,'datadir',[])
    % dataDir = viewGet(view,'datadir')
    if ieNotDefined('varargin')
      groupNum = viewGet(view,'currentGroup');
    else
      groupNum = varargin{1};
    end
    if isempty(groupNum)
      groupNum = viewGet(view,'currentGroup');
    end
    if ~isempty(groupNum)
      viewDir = viewGet(view,'homedir');
      groupname = viewGet(view,'groupname',groupNum);
      val = fullfile(viewDir,groupname);
      if ~exist(val,'dir')
        mkdir(viewDir,groupname);
      end
    end
  case {'roidir'}
    % roiDir = viewGet(view,'roidir')
    viewDir = viewGet(view,'homedir');
    val = fullfile(viewDir,'ROIs');
    if ~exist(val,'dir')
      mkdir(viewDir,'ROIs');
    end
  case {'etcdir'}
    % etcDir = viewGet(view,'etcdir')
    etcDir = fullfile(fileparts(viewGet(view,'datadir',1)),'Etc');
    if isdir(etcDir)
      val = etcDir;
    end
  case {'loadedanalyses'}
    % loadedanalyses = viewGet(view,'loadedAnalyses',[groupNum])
    % this stores all the loaded analyses so that we can switch
    % between groups
    if ieNotDefined('varargin')
      groupNum = viewGet(view,'currentGroup');
    else
      groupNum = varargin{1};
    end
    if isempty(groupNum)
      groupNum = viewGet(view,'currentGroup');
    end
    if ~isempty(groupNum)
      if (groupNum >= 1) & (groupNum <= length(view.loadedAnalyses))
        val = view.loadedAnalyses{groupNum};
      end
    end
  case {'groupscannum'}
    % groupscannum = viewGet(view,'groupscannum',[groupNum])
    % this stores what scan number we were on when we switch
    % between groups
    if ieNotDefined('varargin')
      groupNum = viewGet(view,'currentGroup');
    else
      groupNum = varargin{1};
    end
    if isempty(groupNum)
      groupNum = viewGet(view,'currentGroup');
    end
    if ~isempty(groupNum)
      if (groupNum >= 1) && (groupNum <= length(view.groupScanNum))
        val = view.groupScanNum(groupNum);
      end
    end
    % default to scan 1
    if isempty(val),val = 1;end
  case {'tseriesdir'}
    % tseriesdir = viewGet(view,'tseriesdir',[groupNum])
    % tseriesdir = viewGet(view,'tseriesdir',[])
    % tseriesdir = viewGet(view,'tseriesdir')
    if ieNotDefined('varargin')
      groupNum = viewGet(view,'currentGroup');
    else
      groupNum = varargin{1};
    end
    if isempty(groupNum)
      groupNum = viewGet(view,'currentGroup');
    end
    if ~isempty(groupNum)
      dataDir = viewGet(view,'datadir',groupNum);
      val = fullfile(dataDir,'TSeries');
      if ~exist(val,'dir')
        mkdir(dataDir,'TSeries');
      end
    end
  case {'analysisdir'}
    % analysisdir = viewGet(view,'analysisdir',[groupNum],[analysisNum])
    % analysisdir = viewGet(view,'analysisdir',[groupNum],[])
    % analysisdir = viewGet(view,'analysisdir',[],[analysisNum])
    % analysisdir = viewGet(view,'analysisdir',[],[])
    % analysisdir = viewGet(view,'analysisdir')
    if ieNotDefined('varargin')
      groupNum = viewGet(view,'currentGroup');
      analysisNum = viewGet(view,'currentAnalysis');
    end
    switch length(varargin)
      case 1
        groupNum = varargin{1};
        analysisNum = viewGet(view,'currentAnalysis');
      case 2
        groupNum = varargin{1};
        analysisNum = varargin{2};
    end
    if isempty(groupNum)
      groupNum = viewGet(view,'currentGroup');
    end
    if isempty(analysisNum)
       mrWarnDlg('(viewGet) No analysis is loaded');
       return;
    end
    analysis = view.analyses{analysisNum};
    if ~isempty(groupNum) & ~isempty(analysis)
      dataDir = viewGet(view,'datadir',groupNum);
      val = fullfile(dataDir,analysis.type);
      if ~exist(val,'dir')
        mkdirSuccess = mkdir(dataDir,analysis.type);
	if ~mkdirSuccess
	  mrErrorDlg(sprintf('(viewGet:analysisDir) Could not make directory %s',dataDir));
	end
      end
    end
  case {'overlaydir'}
    % overlaydir = viewGet(view,'overlaydir',[groupNum],[analysisNum])
    val = viewGet(view,'analysisdir',varargin{:});
    % group
  case{'numberofgroups','numgroups','ngroups'}
    % n = viewGet(view,'numberofGroups')
    val = length(MLR.groups);
  case {'groupnames'}
    % groupNames = viewGet(view,'groupNames')
    val = {MLR.groups(:).name};
  case{'currentgroup','curgroup','selectedgroup'}
    % groupNum = viewGet(view,'currentGroup')
    val = view.curGroup;
  case{'group'}
    % groupStructure = viewGet(view,'group',[groupNum])
    % groupStructure = viewGet([],'group',groupNum)
    % groupStructure = viewGet(view,'group',[])
    % groupStructure = viewGet(view,'group')
    if ieNotDefined('varargin')
      groupNum = viewGet(view,'currentGroup');
    else
      groupNum = varargin{1};
    end
    if isempty(groupNum)
      groupNum = viewGet(view,'currentGroup');
    end
    val = MLR.groups(groupNum);
  case{'groupname'}
    % groupName = viewGet(view,'groupName',[groupNum])
    % groupName = viewGet([],'groupName',groupNum)
    % groupName = viewGet(view,'groupName',[])
    % groupName = viewGet(view,'groupName')
    if ieNotDefined('varargin')
      groupNum = viewGet(view,'currentGroup');
    elseif isstr(varargin{1})
      % if passed in a string, look for the string
      % as a group name
      groupNum = viewGet(view,'groupNum',varargin{1});
    else
      groupNum = varargin{1};
    end
    if isempty(groupNum)
      groupNum = viewGet(view,'currentGroup');
    end
    val = MLR.groups(groupNum).name;
  case{'groupnum'}
    % groupnum = viewGet(view,'groupnum',groupName)
    if isstr(varargin{1})
      groupName = varargin{1};
      groupNames = {MLR.groups(:).name};
      val = find(strcmp(groupName,groupNames));
      if isempty(val)
        disp(sprintf('(viewGet) Could not find group: %s',groupName));
      end
      % if passed in a valid number just return that number
    elseif isnumeric(varargin{1}) && isequal(size(varargin{1}),[1 1])
      if (varargin{1} >= 1) && (varargin{1} <= viewGet(view,'nGroups'))
        val = varargin{1};
      end
    end
    % scan
  case{'nscans','numberofscans','numscans'}
    % n = viewGet(view,'nScans',[groupNum])
    % n = viewGet([],'nScans',groupNum)
    % n = viewGet(view,'nScans',[])
    % n = viewGet(view,'nScans')
    g = getGroup(view,varargin);
    if isempty(g)
      val = [];
    else
      val = length(MLR.groups(g).scanParams);
    end
  case{'scanparams'}
    % n = viewGet(view,'scanParams',scanNum,[groupNum])
    % n = viewGet([],'scanParams',scanNum,groupNum)
    % n = viewGet(view,'scanParams',scanNum,[])
    % n = viewGet(view,'scanParams',scanNum)
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      val = MLR.groups(g).scanParams(s);
    end
  case{'auxparams'}
    % n = viewGet(view,'auxParams',scanNum,[groupNum])
    % get the whole auxParams structure. Compare
    % with auxparam which gets one field of auxParams
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      val = MLR.groups(g).auxParams(s);
    end
  case{'stimfile'}
    % stimfile = viewGet(view,'stimfile',[scanNum],[groupNum],[removeEmpties]);
    % returns a cell array with all the stimfiles associated with the
    % scan (including if it is a concatenation or an average, will show
    % all ths stimfiles from the original scans used to make the scan).
    % If you set removeEmpties (defaults to true) to false, then it
    % will return empties for any original scan that has a missing
    % stimfile
    [s g] = getScanAndGroup(view,varargin,param);
    % get removeEp=mpties
    if length(varargin) >= 3,removeEmpties = varargin{3};else removeEmpties = true;end
    % get the stimFileNames from auxParam, calling viewGetLoadStimFileile
    % to prepend etc directory and load
    val = viewGet(view,'auxParam','stimFileName',s,g,@viewGetLoadStimFile,removeEmpties);
    if ~isempty(val),val = cellArray(val);end
  case{'stimfilename'}
    % stimFileName = viewGet(view,'stimFileName',[scanNum],[groupNum],[removeEmpties]);
    % returns the names of all the stimfiles associated with the scan
    % this will prepend the ETc directory where they should live so
    % gives a fully qualified path. If it is a concatenation or an average,
    % will dispaly all stimfile names from the original scans used to make the scan.
    % If you set removeEmpties (defaults to true) to false, then it
    % will return empties for any original scan that has a missing
    % stimfile
    [s g] = getScanAndGroup(view,varargin,param);
    % get the stimfilenames from auxParam, calling viewGetPrependEtc
    % to prepend etc directory
    val = viewGet(view,'auxParam','stimFileName',s,g,@viewGetPrependEtc,1);
    if ~isempty(val),val = cellArray(val);else val = {};end
  case{'fidinfo'}
    % fidinfo = viewGet(view,'fidInfo',[scanNum],[groupNum],[removeEmpties]);
    % returns a cell array with the fidInfo associated with each scan
    % this is only for varian systems in which there is a fid file
    % linked with the auxParam 'fidFilename'
    % (if it is a concatenation or an average, will show
    % all ths fidInfo from the original scans used to make the scan).
    % If you set removeEmpties (defaults to false) to true, then it
    % will not return empties for any original scan that has a missing
    % fidFilename or missing fid
    [s g] = getScanAndGroup(view,varargin,param);
    % get removeEp=mpties
    if length(varargin) >= 3,removeEmpties = varargin{3};else removeEmpties = false;end
    % get the fidFilename from auxParam, calling viewGetFidInfo
    % to prepend Pre directory and get fidInfo
    val = viewGet(view,'auxParam','fidFilename',s,g,@viewGetFidInfo,removeEmpties);
    if ~isempty(val),val = cellArray(val);end
  case{'fidfilename'}
    % fidFileName = viewGet(view,'fidFileName',[scanNum],[groupNum],[removeEmpties]);
    % returns the names of all the fidfiles associated with the scan
    % this will prepend the Pre directory where they should live so
    % gives a fully qualified path. If it is a concatenation or an average,
    % will dispaly all fidfile names from the original scans used to make the scan.
    % If you set removeEmpties (defaults to false) to true, then it
    % will not return empties for any original scan that has a missing
    % fidFileName
    [s g] = getScanAndGroup(view,varargin,param);
    % get the fidFilenames from auxParam, calling viewGetPrependPre
    % to prepend Pre directory
    val = viewGet(view,'auxParam','fidFilename',s,g,@viewGetPrependPre,1);
    if ~isempty(val),val = cellArray(val);else val = {};end
  case{'auxparam'}
    % n = viewGet(view,'auxParam','paramName',[scanNum],[groupNum],[functionptr],<removeEmpties>)
    % get the named parameter from auxParams. Compare
    % with auxParams which gets the whole structure. Note that
    % this will return a single value if paramName exists for that
    % scan. If it doesn't and exist in the original scan (i.e.
    % for a concat or average), then it will return a cell array
    % of the parameter - one for each of the original scans. <functionptr> is
    % an argument used by other viewGet commands - it calls the function
    % on the auxParam field before returning it.
    if length(varargin) < 1
      disp(sprintf('(viewGet) Need to specify paramName for auxParam'));
      return
    end
    [s g] = getScanAndGroup(view,{varargin{2:end}});
    if isempty(g),g = viewGet(view,'curGroup');end
    nscans = viewGet(view,'nscans',g);
    val = [];
    paramName = fixBadChars(varargin{1});
    % handle auxillary arguments
    if length(varargin) >= 4,fhandle = varargin{4};else fhandle = [];end
    if length(varargin) >= 5,removeEmpties = varargin{5};else removeEmpties = false;end
    if (nscans >= s) && (s > 0) && isfield(MLR.groups(g),'auxParams') && (length(MLR.groups(g).auxParams) >= s)
      if isfield(MLR.groups(g).auxParams(s),paramName)
        if ~isempty(MLR.groups(g).auxParams(s).(paramName))
          % get the stored stimFileNames
          val = MLR.groups(g).auxParams(s).(paramName);
	  % val can still be a cell array at this point - usually
	  % because the scan has been copied from somewhere else
	  % and it doesn't have the original raw scan to go back to.
	  % in this case val is a cell array
	  if iscell(val) && (length(val) == 1)
	    val = val{1};
	  end
	  % call function on value (if function handle is specified)
	  if ~isempty(fhandle) 
	    if iscell(val)
	      for i = 1:length(val)
		val{i} = feval(fhandle,view,val{i}); 
	      end
	    else
	      val = feval(fhandle,view,val); 
	    end
	  end
	end
      end
    end
    % if the field does not exist, then check original
    if isempty(val)
      [os og] = viewGet(view,'originalScanNum',s,g);
      if ~isempty(os)
        for osnum = 1:length(os)
          val{osnum} = viewGet(view,'auxParam',paramName,os(osnum),og(osnum),fhandle,removeEmpties);
	  % don't let the cell arrays get too nested
	  if iscell(val{osnum}) && (length(val{osnum}) == 1)
	    val{osnum} = val{osnum}{1};
	  end
        end
	% remove empties
	if removeEmpties
	  newval = {};
	  for j = 1:length(val)
	    if ~isempty(val{j})
	      newval{end+1} = val{j};
	    end
	  end
	  val = newval;
	end
      end
    end
  case{'totalframes'}
    % n = viewGet(view,'totalFrames',scanNum,[groupNum])
    % n = viewGet([],'totalFrames',scanNum,groupNum)
    % n = viewGet(view,'totalFrames',scanNum,[])
    % n = viewGet(view,'totalFrames',scanNum)
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      val = MLR.groups(g).scanParams(s).totalFrames;
    end
  case{'junkframestotal','totaljunkedframes'}
    % gets all junk frames from this scan and the scans this was made from
    % that is the total of all frames that have been junked. This
    % does not included the number of frames in the junkFrames parameter
    % n = viewGet(view,'totalJunkedFrames',scanNum,[groupNum])
    % n = viewGet([],'totalJunkedFrames',scanNum,groupNum)
    % n = viewGet(view,'totalJunkedFrames',scanNum,[])
    % n = viewGet(view,'totalJunkedFrames',scanNum)
    [s g] = getScanAndGroup(view,varargin,param);
    % get the frames that _have already been_ junked from this scan
    val = [];
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0) & isfield(MLR.groups(g).scanParams(s),'totalJunkedFrames')
      val = MLR.groups(g).scanParams(s).totalJunkedFrames;
    end
    if isempty(val)
      % figure out how many time series
      numTimeSeries = max(1,length(viewGet(view,'originalFilename',s,g)));
      % add just set the totalJunkedFrames to 0 for each one
      val = zeros(1,numTimeSeries);
    end
  case{'junkframes'}
    % n = viewGet(view,'junkFrames',scanNum,[groupNum])
    % n = viewGet([],'junkFrames',scanNum,groupNum)
    % n = viewGet(view,'junkFrames',scanNum,[])
    % n = viewGet(view,'junkFrames',scanNum)
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      val = MLR.groups(g).scanParams(s).junkFrames;
    end
  case{'nframes','numframes','nvolumes','numvolumes'}
    % n = viewGet(view,'nFrames',scanNum,[groupNum])
    % n = viewGet([],'nFrames',scanNum,groupNum)
    % n = viewGet(view,'nFrames',scanNum,[])
    % n = viewGet(view,'nFrames',scanNum)
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      val = MLR.groups(g).scanParams(s).nFrames;
    end
  case{'frameperiod'}
    % Frame period is the time from one volume to the next.
    % It is initially extracted from the nifti header. It
    % is *not* necessarily the TR.
    % framePeriod = viewGet(view,'framePeriod',scanNum,[groupNum])
    % framePeriod = viewGet([],'framePeriod',scanNum,groupNum)
    % framePeriod = viewGet(view,'framePeriod',scanNum,[])
    % framePeriod = viewGet(view,'framePeriod',scanNum)
    [s g] = getScanAndGroup(view,varargin,param);
    if isempty(g)
      g = viewGet(view,'currentGroup');
    end
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      val = MLR.groups(g).scanParams(s).framePeriod;
    end
  case{'tseriespathstr','tseriespath'}
    % tseriesPath = viewGet(view,'tseriesPath',scanNum,[groupNum])
    % tseriesPath = viewGet([],'tseriesPath',scanNum,groupNum)
    % tseriesPath = viewGet(view,'tseriesPath',scanNum,[])
    % tseriesPath = viewGet(view,'tseriesPath',scanNum)
    [s g] = getScanAndGroup(view,varargin,param);
    ngroups = viewGet(view,'ngroups');
    if (g > 0) & (g <= ngroups)
      nscans = viewGet(view,'nscans',g);
      if (nscans >= s) & (s > 0)
	val = fullfile(viewGet(view,'tseriesDir',g), MLR.groups(g).scanParams(s).fileName);
      end
    end
  case{'tseriesfile'}
    % filename = viewGet(view,'tseriesFile',scanNum,groupNum)
    % filename = viewGet([],'tseriesFile',scanNum,groupNum)
    [s g] = getScanAndGroup(view,varargin,param);
    if s <= length(MLR.groups(g).scanParams)
      val = MLR.groups(g).scanParams(s).fileName;
    else
      val = [];
    end
  case {'scannum'}
    % get scanNum (and groupNum) of a scan by its file name
    % [scanNum,groupNum] = viewGet(view,'scanNum',tseriesFileName);
    % returns empty if the filename does not exist
    %
    % if you pass in a groupNum as well, then search for
    % tfilename is only done w/in the group.
    % [scanNum,groupNum] = viewGet(view,'scanNum',tseriesFileName,groupNum);
    if ieNotDefined('varargin')
      mrErrorDlg('(viewGet) scanNum: Must specify tseriesFileName.');
    end
    [tseriesFilePath,tseriesFileName] = fileparts(varargin{1});
    if length(varargin) > 1
      groupNums = viewGet(view,'groupNum',varargin{2});
      if isempty(groupNums)
        disp(sprintf('(viewGet) scanNum: Unknown group'));
        return
      end
    else
      groupNums = 1:length(MLR.groups);
    end
    scanNumMatch = [];groupNumMatch = [];
    for groupNum = groupNums
      for scanNum = 1:length(MLR.groups(groupNum).scanParams)
        [scanFilePath,scanFileName] = fileparts(MLR.groups(groupNum).scanParams(scanNum).fileName);
        if strcmp(scanFileName,tseriesFileName)
          scanNumMatch(end+1) = scanNum;
          groupNumMatch(end+1) = groupNum;
        end
      end
    end
    val = scanNumMatch;
    val2 = groupNumMatch;
    % check for more than one scan in the same group that
    % has the same filename. This should never happen since
    % filenames are meant to be unique identifiers of the scan
    uniqueGroupNumMatch = unique(groupNumMatch);
    if length(uniqueGroupNumMatch) ~= length(groupNumMatch)
      % init new arrays
      val = [];
      val2 = [];
      % go through all unique group nums
      for i = 1:length(uniqueGroupNumMatch)
	% find which ones are duplicated so we can display
	whichMatch = find(uniqueGroupNumMatch(i) == groupNumMatch);
	if length(whichMatch) > 1
	  % display warning
	  mrWarnDlg(sprintf('(viewGet:scanNum) !!! Scans %s in group %s have the same filename: %s. Returning only scan %i. !!!\nNote this should not happen because tseries names for scans are meant to be unique identifiers. You should probably delete one of the scans.',num2str(scanNumMatch(whichMatch),'%i '),viewGet(view,'groupName',groupNumMatch(whichMatch(1))),viewGet(view,'tSeriesFile',scanNumMatch(whichMatch(1)),groupNumMatch(whichMatch(1))),scanNumMatch(whichMatch(1))));
	end
	% now keep only the first match
	val(end+1) = scanNumMatch(whichMatch(1));
	val2(end+1) = groupNumMatch(whichMatch(1));
      end
    end
    

  case {'scannumfromdescription'}
    % scanNum = viewGet(view,'scanNumFromDescription',description,<groupNum>,<searchType>);
    % get scanNum(s) in current group that have a matching description
    % returns empty if there are no matches.
    % searchType can be any one of:
    %   'anywhere': case insensitive match where the description string can be found anywhere in the
    %      scans description. e.g. a scan with description 'CCW wedges' will match 'ccw'. Be aware
    %      that 'CCW wedges' will also match the string 'CW'.
    %   'anywhereCaseSensitive': like above, but case sensitive.
    %   'exact': case sensitive match where the description string has to exactly match
    %      the scans description. e.g. 'CCW wedges' will only match 'CCW wedges', not 'CCW'
    %   'exactCaseInsensitive': like exact, but not case sensitive.
    %   'regexp': Will use matlab's regular expression function (regexp) to do the match. For
    %      instance, if you want to match 'CW wedges', but not 'CCW wedges', you could
    %      search for '^CW '.
    % The default serachType is 'anywhere'
    if ieNotDefined('varargin')
      disp('(viewGet) scanNumFromDescription: Must specify description.');
      return
    end
    % get the search string
    searchString = varargin{1};
    % get the group num
    if length(varargin)>1
      groupNum = viewGet(view,'groupNum',varargin{2});
    else
      groupNum = viewGet(view,'curGroup');
    end
    % get the serach Type
    if length(varargin)>2
      searchType = varargin{3};
    else
      searchType = 'anywhere';
    end
    if ~any(strcmp(searchType,{'anywhere','anywhereCaseSensitive','exact','exactCaseInsensitive','regexp'}))
      disp(sprintf('(viewGet) scanNumFromDescription: Unknwon searchType %s',searchType));
      return
    end
    % now go do the search
    for scanNum = 1:viewGet(view,'nScans',groupNum)
      thisDescription = viewGet(view,'description',scanNum,groupNum);
      switch searchType
        case 'anywhere'
          if ~isempty(strfind(lower(thisDescription),lower(searchString)));
            val(end+1) = scanNum;
          end
        case 'anywhereCaseSensitive'
          if ~isempty(strfind(thisDescription,searchString))
            val(end+1) = scanNum;
          end
        case 'exact'
          if strcmp(thisDescription,searchString)
            val(end+1) = scanNum;
          end
        case 'exactCaseInsensitive'
          if strcmp(lower(thisDescription),lower(searchString))
            val(end+1) = scanNum;
          end
        case 'regexp'
          if ~isempty(regexp(thisDescription,searchString))
            val(end+1) = scanNum;
          end
      end
    end
  case {'concatinfo'}
    % concatInfo = viewGet(view,'concatInfo',[scanNum],[groupNum]);
    [s g] = getScanAndGroup(view,varargin,param);
    [tseriesPath,tseriesFile] = fileparts(viewGet(view,'tseriesPath',s,g));
    % check for mat file
    matFileName = fullfile(tseriesPath,sprintf('%s.mat',tseriesFile));
    if isfile(matFileName)
      load(matFileName);
      if exist('concatInfo')
        val = concatInfo;
      end
    end
  case {'stimfileold'}
    % stimfile = viewGet(view,'stimfile',scanNum,[groupNum]);
    % this is the old way of getting stimfiles. The new version
    % should be completely compatible with this old version,
    % but leaving this code in here, just in case, for now.
    % came be deleted later. jg 11/2011
    [s g] = getScanAndGroup(view,varargin,param);
    stimFileName = viewGet(view,'stimFileName',s,g);
    val = {};
    % cycle over all stimfiles (concatenations and averages,
    % may have several stimfiles).
    if ~isempty(stimFileName)
      for j = 1:length(stimFileName)
        % load this stimfile
        if ~isfile(stimFileName{j})
          mrErrorDlg(sprintf('viewGet %s: Could not find stimfile %s',param,stimFileName{j}));
        else
          val{j}=load(stimFileName{j});
        end
        % check to see what type it is, and set the field appropriately
        if isfield(val{j},'mylog')
          val{j}.filetype = 'eventtimes';
        elseif isfield(val{j},'stimts')
          val{j}.filetype = 'afni';
        elseif isfield(val{j},'myscreen')
          val{j}.filetype = 'mgl';
        elseif isfield(val{j},'stimvol')
          val{j}.filetype = 'stimvol';
        else
          val{j}.filetype = 'unknown';
        end
      end
    end
  case {'stimfilenameold'}
    % stimFileName = viewGet(view,'stimFileName',scanNum,[groupNum]);
    % stimFileName is returned as a cell array of all stimfiles
    % associated with the scan
    % this is the old way of getting stimfiles. The new version
    % should be completely compatible with this old version,
    % but leaving this code in here, just in case, for now.
    % came be deleted later. jg 11/2011
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    stimFileName{1} = '';
    if (nscans >= s) && (s > 0) && isfield(MLR.groups(g),'auxParams') && (length(MLR.groups(g).auxParams) >= s)
      if isfield(MLR.groups(g).auxParams(s),'stimFileName')
        if ~isempty(MLR.groups(g).auxParams(s).stimFileName)
          % get the stored stimFileNames
          stimFileName = cellArray(MLR.groups(g).auxParams(s).stimFileName);
          % and prepend the etcDir path on to them
          for i = 1:length(stimFileName)
            stimFileName{i} = fullfile(viewGet(view,'etcDir'),stimFileName{i});
          end
        end
      end
    end
    % if the file does not exist, then check original
    if isempty(stimFileName{1})
      stimFileName = {};
      [os og] = viewGet(view,'originalScanNum',s,g);
      if ~isempty(os)
        for osnum = 1:length(os)
          newStimFileNames = viewGet(view,'stimFileName',os(osnum),og(osnum));
          for j = 1:length(newStimFileNames)
            stimFileName{end+1} = newStimFileNames{j};
          end
        end
      end
    end
    val = stimFileName;
  case{'mousedownscancoords'}
    % scanCoords = viewGet(view,'mouseDownScanCoords')
    viewNum = viewGet(view,'viewNum');
    if isfield(MLR,'interrogator')
      if length(MLR.interrogator) >= viewNum
        if isfield(MLR.interrogator{viewNum},'mouseDownScanCoords')
          val = MLR.interrogator{viewNum}.mouseDownScanCoords;
        end
      end
    end
  case{'mousedownbasecoords'}
    % baseCoords = viewGet(view,'mouseDownBaseCoords')
    viewNum = viewGet(view,'viewNum');
    if isfield(MLR,'interrogator')
      if length(MLR.interrogator) >= viewNum
        if isfield(MLR.interrogator{viewNum},'mouseDownBaseCoords')
          val = MLR.interrogator{viewNum}.mouseDownBaseCoords;
        end
      end
    end
  case {'spikeinfo'}
    % eyepos = viewGet(view,'spikeinfo',scanNum,[groupNum]);
    val = [];
    [s g] = getScanAndGroup(view,varargin,param);
    nScans = viewGet(view,'nscans',g);
    nGroups = viewGet(view,'nGroups');
    if (g > 0) &  (length(MLR.groups(g).auxParams) >= s) & (s > 0)
      if isfield(MLR.groups(g).auxParams(s),'spikeinfo')
        val = MLR.groups(g).auxParams(s).spikeinfo;
        % check tSeries name
        if isfield(val,'filename') & ~strcmp(val.filename,viewGet(view,'tSeriesFile',s,g))
          val = [];
        end
      end
    end
  case {'eyepos'}
    % eyepos = viewGet(view,'eyepos',scanNum,[groupNum]);
    [s g] = getScanAndGroup(view,varargin,param);
    [eyeposFileName eyeposNum] = viewGet(view,'eyeposFileName',s,g);
    val = {};
    % cycle over all eyepos files (concatenations and averages,
    if ~isempty(eyeposFileName)
      for j = 1:length(eyeposFileName)
        % load this eyepos file
        if ~isfile(sprintf('%s.mat',stripext(eyeposFileName{j})))
          mrErrorDlg(sprintf('viewGet %s: Could not find eyepos file %s',param,eyeposFileName{j}));
        else
          eyepos=load(eyeposFileName{j});
          f = fieldnames(eyepos);
          % here we get the proper scan.
          val{j} = eyepos.(f{1}).scan{eyeposNum{j}};
          val{j}.hdr = rmfield(eyepos.(f{1}),'scan');
        end
      end
    end
  case {'eyeposfilename'}
    % eyeposFileName = viewGet(view,'eyeposFileName',scanNum,[groupNum]);
    % eyeposFileName is returned as a cell array of all eyepos files
    % associated with the scan
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    eyeposFileName{1} = '';
    if (nscans >= s) && (s > 0) && isfield(MLR.groups(g),'auxParams') && (length(MLR.groups(g).auxParams) >= s)
      if isfield(MLR.groups(g).auxParams(s),'eyeposFileName')
        if ~isempty(MLR.groups(g).auxParams(s).eyeposFileName)
          eyeposFileName{1} = fullfile(viewGet(view,'etcDir'),MLR.groups(g).auxParams(s).eyeposFileName);
          eyeposNum{1} = MLR.groups(g).auxParams(s).eyeposNum;
        end
      end
    end
    % if the file does not exist, then check original
    if isempty(eyeposFileName{1})
      eyeposFileName = {};
      eyeposNum = {};
      [os og] = viewGet(view,'originalScanNum',s,g);
      if ~isempty(os)
        for osnum = 1:length(os)
          [newEyeposFileNames newEyeposNum] = viewGet(view,'eyeposFileName',os(osnum),og(osnum));
          for j = 1:length(newEyeposFileNames)
            eyeposFileName{end+1} = newEyeposFileNames{j};
            eyeposNum{end+1} = newEyeposNum{j};
          end
        end
      end
    end
    val = eyeposFileName;
    val2 = eyeposNum;
  case {'3d'}
    % is3d = viewGet(view,'tr',scanNum,[groupNum]);
    % returns whether sequence is 3d or not
    [s g] = getScanAndGroup(view,varargin,param);
    % first get dicom
    dicom = viewGet(view,'dicom',s,g);
    val = [];
    if isempty(dicom)
      return
    end
    for i = 1:length(dicom)
      if isfield(dicom{i},'MRAcquisitionType')
	thisval = dicom{i}.MRAcquisitionType;
      elseif isfield(dicom{i},'ACQ') && isfield(dicom{i}.ACQ,'MR_Acquisition_Type_')
        thisval = strcmp(dicom{i}.ACQ.MR_Acquisition_Type_,'3D');
      else
	thisval = [];
      end
      if (isempty(val))
        val = thisval;
      elseif (val ~= thisval)
        disp(sprintf('(viewGet:3D) Data collected with both 3D and 2D sequence types'))
      end
    end
  case {'tr'}
    % TR is the TR of the pulse sequence. It is *not* necessarily
    % the time between volumes (framePeriod). The TR is extracted
    % from the DICOM header
    % tr = viewGet(view,'tr',scanNum,[groupNum]);
    [s g] = getScanAndGroup(view,varargin,param);
    % first get dicom
    dicom = viewGet(view,'dicom',s,g);
    val = [];
    if isempty(dicom)
      return
    end
    for i = 1:length(dicom)
      thistr = [];
      if isfield(dicom{i},'RepetitionTime')
	thistr = dicom{i}.RepetitionTime/1000;
      elseif isfield(dicom{i},'ACQ') && isfield(dicom{i}.ACQ,'Repetition_Time')
        thistr = dicom{i}.ACQ.Repetition_Time/1000;
      end
      if  ~isempty(thistr)
        if (isempty(val))
          val = thistr;
        elseif (val ~= thistr)
          disp(sprintf('(viewGet) Data collected with different TR: %0.2f vs %0.2f',val,thistr));
        end
      end
    end
  case {'dicom'}
    % dicom = viewGet(view,'dicom',scanNum,[groupNum]);
    [s g] = getScanAndGroup(view,varargin,param);
    % see if we have dicom info in auxparam
    val = viewGet(view,'auxParam','dicomInfo',s,g);
    if isempty(val)
      dicomName = viewGet(view,'dicomName',s,g);
      val = {};
      if ~isempty(dicomName)
	for j = 1:length(dicomName)
	  val{j} = readdicomheader(dicomName{j});
	end
      end
    end
    % make into a cell array
    if ~isempty(val)
      val = cellArray(val);
    end
  case {'dicomname'}
    % dicomName = viewGet(view,'dicomName',scanNum,[groupNum]);
    % dicomName is returned as a cell array of all dicom headers
    % associated with the scan
    [s g] = getScanAndGroup(view,varargin,param);
    % create name of file, by putting -header.txt on to the filename
    [tseriesPath,tseriesFile] = fileparts(viewGet(view,'tseriesPath',s,g));
    dicomName{1} = fullfile(tseriesPath,sprintf('%s-header.txt',tseriesFile));
    % if the file does not exist, then check original
    if ~isfile(dicomName{1})
      dicomName = {};
      [os og] = viewGet(view,'originalScanNum',s,g);
      if ~isempty(os)
        for osnum = 1:length(os)
          newDicomNames = viewGet(view,'dicomName',os(osnum),og(osnum));
          for j = 1:length(newDicomNames)
            dicomName{end+1} = newDicomNames{j};
          end
        end
      end
    end
    val = dicomName;
  case {'originalscannum'}
    % [scanNum,groupNum] = viewGet(view,'originalScanNum',scanNum,[groupNum]);
    [s g] = getScanAndGroup(view,varargin,param);
    originalFileNames = viewGet(view,'originalFileName',s,g);
    originalGroupNames = viewGet(view,'originalGroupName',s,g);
    originalScanNum = [];originalGroupNum = [];
    for i = 1:length(originalFileNames)
      % make sure there is a valid scan to return
      if ~isempty(viewGet(view,'scanNum',originalFileNames{i},viewGet(view,'groupNum',originalGroupNames{i})))
        [originalScanNum(end+1),originalGroupNum(end+1)] = viewGet(view,'scanNum',originalFileNames{i},viewGet(view,'groupNum',originalGroupNames{i}));
      end
    end
    val = originalScanNum;val2 = originalGroupNum;
  case {'originalfilename'}
    % originalFileName = viewGet(view,'originalFileName',scanNum,[groupNum]);
    % if the original is the same as the filename, then this returns []
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      if isfield(MLR.groups(g).scanParams(s),'originalFileName')
        val = MLR.groups(g).scanParams(s).originalFileName;
        % if originalFilename is the same as the filename
        % then return empty
        if strcmp(val,MLR.groups(g).scanParams(s).fileName)
          val = [];
        end
      else
        val = [];
      end
    end
  case {'originalgroupname'}
    % originalGroupName = viewGet(view,'originalGroupName',scanNum,[groupNum]);
    % if the original is the same as the filename, then this returns []
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      if isfield(MLR.groups(g).scanParams(s),'originalGroupName')
        val = MLR.groups(g).scanParams(s).originalFileName;
        if strcmp(val,MLR.groups(g).scanParams(s).fileName)
          val = [];
        else
          val = MLR.groups(g).scanParams(s).originalGroupName;
        end
      else
        val = [];
      end
    end
  case{'transforms'}
    % transforms = viewGet(view,'transforms',scanNum,[groupNum]);
    % returns motion correction transformation matrices if they exists
    [s g] = getScanAndGroup(view,varargin,param);
    % get the tseries name
    [tSeriesPath tSeriesName] = fileparts(viewGet(view,'tSeriesPathStr',s,g));
    if isfile(fullfile(tSeriesPath,sprintf('%s.mat',tSeriesName)))
      load(fullfile(tSeriesPath,sprintf('%s.mat',tSeriesName)));
      if exist('transforms') == 1
        val = transforms;
      end
    end
  case{'params'}
    % params = viewGet(view,'params',scanNum,[groupNum]);
    % gets the .mat params file associated with this scan if it exists
    [s g] = getScanAndGroup(view,varargin,param);
    % get the tseries name
    [tSeriesPath tSeriesName] = fileparts(viewGet(view,'tSeriesPathStr',s,g));
    if isfile(fullfile(tSeriesPath,sprintf('%s.mat',tSeriesName)))
      val = load(fullfile(tSeriesPath,sprintf('%s.mat',tSeriesName)));
    end
  case{'niftihdr'}
    % hdr = viewGet(view,'niftiHdr',scanNum,[groupNum])
    % hdr = viewGet([],'niftiHdr',scanNum,groupNum)
    % hdr = viewGet(view,'niftiHdr',scanNum,[])
    % hdr = viewGet(view,'niftiHdr',scanNum)
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      % make sure it is a nifti headr (not an analyze)
      val = MLR.groups(g).scanParams(s).niftiHdr;
    end
  case{'sformcode','scansformcode'}
    % xform = viewGet(view,'sformcode',[scanNum],[groupNum])
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      if ~isfield(MLR.groups(g).scanParams(s).niftiHdr,'sform_code')
        val = [];
      else
        val = MLR.groups(g).scanParams(s).niftiHdr.sform_code;
      end
    end
  case{'scanqformcode'}
    % xform = viewGet(view,'qformCode',[scanNum],[groupNum])
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      if ~isfield(MLR.groups(g).scanParams(s).niftiHdr,'qform_code')
        val = [];
      else
        val = MLR.groups(g).scanParams(s).niftiHdr.qform_code;
      end
    end
  case{'scanvol2tal'}
    % xform = viewGet(view,'scanVol2tal',[scanNum],[groupNum])
    % This will return the xform matrix that specifies the
    % transformation from volume coordinates to talairach
    % coordinates of the base volume that this scan was aligned to.
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      val = MLR.groups(g).scanParams(s).vol2tal;
    end
  case{'scanvol2mag'}
    % xform = viewGet(view,'scanVol2mag',[scanNum],[groupNum])
    % This will return the xform matrix that specifies the
    % transformation from volume coordinates to the magnet
    % coordinates of the base volume that this scan was aligned to.
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      val = MLR.groups(g).scanParams(s).vol2mag;
    end
  case{'scan2tal'}
    % xform = viewGet(view,'scan2tal',[scanNum],[groupNum])
    % This will return the xform matrix that specifies the
    % transformation from this scan to the volume
    % in talairach coordinates. Note that this also has
    % composited the shiftOriginXform.
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      sform_code = viewGet(view,'sformCode',s,g);
      scanSform = viewGet(view,'scanSform',s,g);
      if ~isempty(scanSform)
        if (sform_code == 3)
          val = scanSform * shiftOriginXform;
        elseif (sform_code == 1)
          vol2tal = viewGet(view,'scanVol2tal',s,g);
          vol2mag = viewGet(view,'scanVol2mag',s,g);
          if ~isempty(vol2tal) && ~isempty(vol2mag)
            val = vol2tal * inv(vol2mag) * scanSform * shiftOriginXform;
          end
        end
      end
    end
  case{'scan2mag','scanxform'}
    % xform = viewGet(view,'scan2mag',[scanNum],[groupNum])
    % The scan2mag xform specifies the xform from the scan
    % to the volume in magnet coordinates.
    % If the sform_code is set to 1, then this is the same
    % as scanSform *except* that the origin has been shifted
    % to start at 1,1,1. i.e. It is scanSform * shiftOriginXform
    %
    % Note that before adding the talairach xform code, scan2mag
    % was called scanxform and the code assumed that the sform
    % of the scan contained this xform. If you ask for scanXform
    % you will get the scan2mag xform, but it does *not* have
    % the shiftOriginXform composited so as to be compatible
    % with the old code
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      sform_code = viewGet(view,'sformCode',s,g);
      scanSform = viewGet(view,'scanSform',s,g);
      if ~isempty(scanSform)
        if (sform_code == 1)
          val = scanSform * shiftOriginXform;
        elseif (sform_code == 3)
          vol2tal = viewGet(view,'scanVol2tal',s,g);
          vol2mag = viewGet(view,'scanVol2mag',s,g);
          if ~isempty(vol2tal) && ~isempty(vol2mag)
            val = vol2mag * inv(vol2tal) * scanSform * shiftOriginXform;
          end
        elseif (sform_code == 0)
          % If sform has not been set, then use the transform that
          % transforms this image directly on to the current anatomy
          % using the qform matrices.
          if strcmp(mrGetPref('verbose'),'Yes')
            oneTimeWarning(sprintf('noScanSform_%i_%i',s,g),...
              ['(viewGet:scanXform) sform is not set. Using qform to align '...
              'to base anatomy. Run mrAlign then mrUpdateNiftiHdr to fix this']);
          end
          val = MLR.groups(g).scanParams(s).niftiHdr.qform44 * shiftOriginXform;
        end
      end
      if strcmp(lower(param),'scanxform') && ~isempty(val)
        val = val * inv(shiftOriginXform);
      end
    end
  case{'scan2roi'}
    % xform = viewGet(view,'scan2roi',[roiNum],[scanNum],[groupNum])
    % This will return the xform matrix that specifies the
    % transformation from scan coordinates to ROI coordinates
    % It checks whether the ROI and the scan are in magnet or Talairach
    % coordinates, and deals with the case when they are not in the
    % same space
    %  Note that this also has composited the shiftOriginXform.
    [s g] = getScanAndGroup(view,varargin,param,2);
    nscans = viewGet(view,'nscans',g);
    r = getRoiNum(view,varargin);
    nRois = viewGet(view,'numrois');
    if (nscans >= s) && (s > 0) && (nRois >= r) && (r > 0)
      roi2tal = viewGet(view,'roi2tal',r);
      if ~isempty(roi2tal) % The roi has a Tal xform
        scan2tal = viewGet(view,'scan2tal',s,g); % check scan
        if ~isempty(scan2tal) % -CASE 1-: both the roi and the scan have a Tal xform
          val = inv(roi2tal) * scan2tal; % use it
        else %  -CASE 2-: the roi has a Tal xform but the scan does not
          roi2mag = viewGet(view,'roi2mag',r); % check if the roi has a Mag xform
          if ~isempty(roi2mag) % if it does,
            scan2mag = viewGet(view,'scan2mag',s,g); % check if the scan has a mag xform too
            if ~isempty(scan2mag) % if they both do, use that, but give a warning
              oneTimeWarning(sprintf('roiScanMismatch_%i_%i_%i',r,s,g),...
                ['WARNING: Ignoring the roi talairach transformation, '...
                'because roi has a Talairach transformation, but '...
                'the scan does not. Coordinates are being converted '...
                'through the magnet coordinate frame. If you wanted '...
                'the transformations to use Talairach coordinates '...
                'instead of magnet coordinates, you need to use '...
                'mrAlign to export the talairach transformation to the scan']);
              val = inv(roi2mag) * scan2mag;
            else % if the scan doesn't have either xform, that's an error.
              % Give a warning and use the identity matrix
              oneTimeWarning(sprintf('noScanXform_%i_%i_%i',r,s,g),...
                ['ERROR: Scan does not have a transform for '...
                'magnet or talairach space. '...
                'Using the identity matrix to transform from '...
                'scan to ROI. Run mrAlign to fix the scan.']);
              val = eye(4);
            end
          else % if the roi does not have a mag xform and the scan does not have a tal xform
            oneTimeWarning(sprintf('incompatibleRoiScan_%i_%i_%i',r,s,g),...
              ['Scan and ROI are not compatible: ROI is in Talairach space '...
              '(and not magnet space) but Scan is not. '...
              'Using the identity matrix to transform from scan to ROI. '...
              'Run mrAlign to get scan and ROI into the same space.']);
            val = eye(4);
          end
        end
      else % The ROI doesn't have a Tal xform...
        roi2mag = viewGet(view,'roi2mag',r);
        if ~isempty(roi2mag) % ... but the ROI does have a mag transform
          scan2mag = viewGet(view,'scan2mag',s,g); % check the scan:
          if ~isempty(scan2mag) % -CASE 3-: both scan and ROI have mag transform
            val = inv(roi2mag) * scan2mag; % use it;
            % but check if scan had a Tal xform so can warn user that it's being ignored
            scan2tal = viewGet(view,'scan2tal',s,g);
            if ~isempty(scan2tal)
              oneTimeWarning(sprintf('roiScanMismatch_%i_%i_%i',r,s,g),...
                ['WARNING: Ignoring the scan talairach transformation, '...
                'because the scan has a talairach '...
                'transformation but the ROI does not. Coordinates '...
                'are being converted through the magnet '...
                'coordinate frame. If you want convert using the '...
                'talairach transformations, you need '...
                'to export a talairach transformation to the '...
                'ROI by running  mrAlign.']);
            end
          else % -CASE 4-: ROI has a mag xform but scan does not
            % (and ROI doesn't have a tal xform, bc already checked that)
            oneTimeWarning(sprintf('incompatibleRoiScan_%i_%i_%i',r,s,g),...
              ['ROI and Scan are not compatible: ROI is in magnet '...
              'space and Scan is not. Using the '...
              'identity matrix to transform from scan to ROI '...
              'Run mrAlign to get scan and ROI into the same space.']);
            val = eye(4);
          end
        else % error if ROI has neither a magnet nor a tal transform
          oneTimeWarning(sprintf('unknownSformCode_%i_%i_%i',r,s,g),...
            ['ROI is neither in Magnet space nor in Talairach '...
            'Space. Using the identity matrix '...
            'to transform from base to ROI. Run mrAlign to fix the ROI.']);
          val = eye(4);
        end
      end
    end
  case{'scan2scan'}
    % xform = viewGet(view,'scan2scan',[fromScanNum],[fromGroupNum],[toScanNum],[toGroupNum])
    % This will return the xform matrix that specifies the
    % transformation from coordinates of the 'From' scan to coordinates of the 'To' scan
    % (E.g., given x,y,s from the 'From' Scan, multiply by the xform calculated in this
    % call to convert to x',y',s' in the 'To' scan.)
    % It checks whether the scans are each in magnet or Talairach
    % coordinates, and deals with the case when they are not in the
    % same space
    %  Note that this also has composited the shiftOriginXform.
    [sFrom gFrom] = getScanAndGroup(view,varargin,param);
    [sTo gTo] = getScanAndGroup(view,varargin,param,3);
    nscansFrom = viewGet(view,'nscans',gFrom);
    nscansTo = viewGet(view,'nscans',gTo);
    if (nscansFrom >= sFrom) && (sFrom > 0) && (nscansTo >= sTo) && (sTo > 0)
      scan2talFrom = viewGet(view,'scan2tal',sFrom,gFrom); % check if the FROMscan
      % has a tal xform
      if ~isempty(scan2talFrom) % if the FROMscan has a Tal xform
        scan2talTo = viewGet(view,'scan2tal',sTo,gTo); % check TOscan
        if ~isempty(scan2talTo) % -CASE 1-: both scans have a Tal xform, then use it
          % first check if they are the same, then just
          % return identity matrix, % so there's no roundoff
          %  error when compute the inverse
          if isequal(scan2talTo, scan2talFrom)
            val = eye(4);
          else
            val = inv(scan2talTo) * scan2talFrom;
          end
        else %  -CASE 2-: the FROMscan has a Tal xform but the TOscan doesn't
          scan2magFrom = viewGet(view,'scan2mag',sFrom,gFrom); % check if the FROMscan
          % has a Mag xform
          if ~isempty(scan2magFrom) % if it does,
            scan2magTo = viewGet(view,'scan2mag',sTo,gTo); % check if the TOscan
            % has a Mag xform too
            if ~isempty(scan2magTo) % if they both do, use that, but give a warning
              oneTimeWarning(sprintf('scanMismatch_%i_%i_%i_%i',sFrom,gFrom,sTo,gTo),...
                ['WARNING: Ignoring the scan talairach transformation, '...
                'because FROM scan has a Talairach transformation,'...
                ' but TO scan does not. Coordinates are being converted'...
                ' through the magnet coordinate frame.'...
                ' If you wanted the transformations to use Talairach'...
                ' coordinates instead of magnet coordinates, you'...
                ' need to use mrAlign to export the talairach '...
                ' transformation to the TO scan']);
              if isequal(scan2magTo,scan2magFrom) % if they're the same, avoid
                % the roundoff errors of calling inv
                val = eye(4);
              else
                val = inv(scan2magTo) * scan2magFrom;
              end
            else % if the TOscan doesn't have either xform, that's an error.
              % Give a warning and use the identity matrix
              oneTimeWarning(sprintf('noScanXform_%i_%i_%i_%i',sFrom,gFrom,sTo,gTo),...
                ['ERROR: TO scan does not have a transform for magnet '...
                'or talairach space. Using the identity'...
                ' matrix to transform from scan to scan. '...
                'Run mrAlign to fix the TO scan.']);
              val = eye(4);
            end
          else % if the FROM scan does not have a mag xform
            % and the TO scan does not have a tal xform
            oneTimeWarning(sprintf('incompatibleScans_%i_%i_%i_%i',sFrom,gFrom,sTo,gTo),...
              ['Scans are not compatible: FROM scan is in '...
              'Talairach space but TO scan is not. '...
              'Using the identity matrix to transform from scan to scan. '...
              'Run mrAlign to get scans into the same space.']);
            val = eye(4);
          end
        end
      else % The FROM scan doesn't have a Tal xform...
        scan2magFrom = viewGet(view,'scan2mag',sFrom,gFrom);
        if ~isempty(scan2magFrom) % ... but the FROM scan does have a mag transform
          scan2magTo = viewGet(view,'scan2mag',sTo,gTo); % check the TO scan
          if ~isempty(scan2magTo) % -CASE 3-: both scans have mag transform
            if isequal(scan2magTo,scan2magFrom) % if they're the same, avoid
              % the roundoff errors of calling inv
              val = eye(4);
            else
              val = inv(scan2magTo) * scan2magFrom;
            end
            % but check if TO scan had a Tal xform so can warn user that it's being ignored
            scan2talTo = viewGet(view,'scan2tal',sTo,gTo);
            if ~isempty(scan2talTo)
              oneTimeWarning(sprintf('scanMismatch_%i_%i_%i_%i',sFrom,gFrom,sTo,gTo),...
                ['WARNING: Ignoring the TO scans talairach transformation, '...
                'because the TO scan has a talairach transformation but '...
                'the FROM scan does not. Coordinates are being converted '...
                'through the magnet coordinate frame. '...
                'If you want convert using the talairach transformations, '...
                'you need to export a talairach transformation to '...
                'the TO scan by running  mrAlign.']);
            end
          else % -CASE 4-: FROM scan has a mag xform but TO scan does not
            % (and scan doesn't have a tal xform, bc already checked that)
            oneTimeWarning(sprintf('incompatibleScans_%i_%i_%i_%i',sFrom,gFrom,sTo,gTo),...
              ['Scans are not compatible: FROM scan is in magnet space but '...
              'TO scan is not. Using the identity matrix to transform from '...
              'scan to scan. Run mrAlign to get scans into the same space.']);
            val = eye(4);
          end
        else % error if scan has neither a base nor a tal transform
          oneTimeWarning(sprintf('unknownSformCode_%i_%i_%i_%i',sFrom,gFrom,sTo,gTo),...
            ['FROM scan is neither in Magnet space nor in Talairach Space.'...
            ' Using the identity matrix to transform from scan to scan.'...
            ' Run mrAlign to fix the scan.']);
          val = eye(4);
        end
      end
    end
  case{'scanqform'}
    % sform = viewGet(view,'scanQform',[scanNum],[groupNum])
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      % get the field
      if isfield(MLR.groups(g).scanParams(s).niftiHdr,'qform44')
        val = MLR.groups(g).scanParams(s).niftiHdr.qform44;
      end
    end
  case{'scansform'}
    % sform = viewGet(view,'scanSform',[scanNum],[groupNum])
    % The scanSform is the sform set in the nifti header.
    % If the sform_code is set to 1, then this field has
    % been set by mrAlign to be the tansformation from this
    % scan to the volume in magnet coordinates. Note that
    % scanSform does not shift the origin to start at 1,1,1
    % You usually will need to composite this sform with
    % shiftOriginXform to get the scan2mag xform.
    % If the sform_code is set to 3, then this field has
    % been set by mrAlign to be the transformation from
    % this scan to the volume in talairach coordinates.
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      % make sure it is a nifti headr (not an analyze)
      if ~isfield(MLR.groups(g).scanParams(s).niftiHdr,'sform_code')
        return
      end
      sformCode = MLR.groups(g).scanParams(s).niftiHdr.sform_code;
      if sformCode
        % if sform has been set, then use it.
        val = MLR.groups(g).scanParams(s).niftiHdr.sform44;
      else
        % If has not been set, then use the transform that
        % transforms this image directly on to the current anatomy
        % using the qform matrices. Note: There used to be code
        % here that reset val if it was the identity but that was a
        % bug (DJH 1/17/07).
        if strcmp(mrGetPref('verbose'),'Yes')
          oneTimeWarning(sprintf('noScanSform_%i_%i',s,g),...
            ['(viewGet:scanXform) sform is not set. Using qform to align '...
            'to base anatomy. Run mrAlign then mrUpdateNiftiHdr to fix this']);
        end
        val = MLR.groups(g).scanParams(s).niftiHdr.qform44;
      end
    end
  case{'scanvoxelsize'}
    % voxelSize = viewGet(view,'scanVoxelSize',scanNum,[groupNum])
    % voxelSize = viewGet([],'scanVoxelSize',scanNum,groupNum)
    % voxelSize = viewGet(view,'scanVoxelSize',scanNum,[])
    % voxelSize = viewGet(view,'scanVoxelSize',scanNum)
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      val = MLR.groups(g).scanParams(s).voxelSize;
    end
  case{'scandims'}
    % dims = viewGet(view,'scanDims',scanNum,[groupNum])
    % dims = viewGet([],'scanDims',scanNum,groupNum)
    % dims = viewGet(view,'scanDims',scanNum,[])
    % dims = viewGet(view,'scanDims',scanNum)
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      val = MLR.groups(g).scanParams(s).dataSize;
    end
  case{'description'}
    % string = viewGet(view,'description',scanNum,[groupNum])
    % string = viewGet([],'description',scanNum,groupNum)
    % string = viewGet(view,'description',scanNum,[])
    % string = viewGet(view,'description',scanNum)
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      val = MLR.groups(g).scanParams(s).description;
    end
    
    % baseVolume (anatomy)
  case{'numberofbasevolumes','numbase'}
    % n = viewGet(view,'numberofBaseVolumes')
    val = length(view.baseVolumes);
  case {'curbase','currentbase'}
    % baseNum = viewGet(view,'currentBase')
    val = view.curBase;
  case{'basenum'}
    % baseNum = viewGet(view,'baseNum',baseName)
    baseName = varargin{1};
    % if numeric there is nothing to do, just return value
    if isnumeric(baseName)
      val = baseName;
    else
      % otherwise look up the baseNum
      baseNames = {view.baseVolumes(:).name};
      val = find(strcmp(baseName,baseNames));
    end
  case{'basevolume','baseanatomy'}
    % basevolume = viewGet(view,'baseVolume',[baseNum])
    if ieNotDefined('varargin')
      b = viewGet(view,'currentBase');
    else
      b = varargin{1};
    end
    if isempty(b)
      b = viewGet(view,'currentBase');
    end
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      val = view.baseVolumes(b);
    end
  case {'currentbasename','curbasename'}
    % baseName = viewGet(view,'currentBaseName')
    baseNames = viewGet(view,'baseNames');
    curBase = viewGet(view,'curBase');
    if isempty(curBase)
      val = 'NoBase';
    else
      val = baseNames{curBase};
    end
  case {'basenames'}
    % baseNames = viewGet(view,'baseNames')
    if ieNotDefined('varargin')
      baseNum = viewGet(view,'curBase');
    else
      baseNum = varargin{1};
    end
    if isempty(baseNum)
      baseNum = viewGet(view,'curBase');
    end
    if ~isempty(view.baseVolumes)
      val = {view.baseVolumes(:).name};
    end
  case {'basename'}
    % basename = viewGet(view,'basename',[baseNum])
    b = getBaseNum(view,varargin);
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      val = view.baseVolumes(b).name;
    end
  case {'basedata'}
    % basedata = viewGet(view,'basedata',[baseNum])
    b = getBaseNum(view,varargin);
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      val = view.baseVolumes(b).data;
    end
  case {'basehandle'}
    % basedata = viewGet(view,'baseHandle',[baseNum])
    b = getBaseNum(view,varargin);
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      val = view.baseVolumes(b).h;
      % set to empty if handle is stale 
      if ~ishandle(val),val = [];end
    end
  case {'basehdr'}
    % basedata = viewGet(view,'basehdr',[baseNum])
    b = getBaseNum(view,varargin);
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      val = view.baseVolumes(b).hdr;
    end
  case {'base'}
    % basedata = viewGet(view,'base',[baseNum])
    b = getBaseNum(view,varargin);
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      val = view.baseVolumes(b);
      [tf val] = isbase(val);
    end
  case {'basecoordmappath'}
    % basedata = viewGet(view,'baseCoordMapPath',[baseNum],[corticalDepth])
    b = getBaseNum(view,varargin);
    n = viewGet(view,'numberofbasevolumes');
    val = [];
    if b & (b > 0) & (b <= n)
      if isfield(view.baseVolumes(b),'coordMap')
        if ~isempty(view.baseVolumes(b).coordMap)
          val = view.baseVolumes(b).coordMap.path;
	  % check if the path exists
	  if ~isdir(val)
	    % guess where the path is in from the canonical volume directory
	    volumeDirectory = mrGetPref('volumeDirectory');
	    volumeDirectoryList = dir(volumeDirectory);
	    innerCoordsFilename = view.baseVolumes(b).coordMap.innerCoordsFileName;
	    subjectDir = '';
	    % tell user what we are doing
	    disp(sprintf('(viewGet:baseCoordMapPath) Surface directory %s for base %s does not exist, searching in volumeDirectory: %s',val,viewGet(view,'baseName',b),volumeDirectory));
	    for i = 1:length(volumeDirectoryList)
	      % for each volume directory in the list, see if the directory name
	      % matches the first part of the baseVolumes anatomy (this assumes
	      % that people use a convention like calling the directory s001 and
	      % calling the anatomy file s001anatomy or something like that.
	      matchName = strfind(view.baseVolumes(b).coordMap.anatFileName,volumeDirectoryList(i).name);
	      if ~isempty(matchName) && isequal(matchName(1),1)
		% we have a match, for the subject directory under the volume direcotry
		subjectDir = fullfile(volumeDirectory,volumeDirectoryList(i).name);
		break;
	      end
	    end
	    % not found, give up
	    if isempty(subjectDir)
	      disp(sprintf('(viewGet:baseCoordMapPath) Could not find a matchind subjectDir in %s',volumeDirectory));
	    else
	      % set to return subjectDir in case we don't find the surfRelax directory
	      val = subjectDir;
	      % Check for a directory right under this one called any of the following
	      possibleSurfDirNames = {'surfRelax','surfaces','surface','Surface','Surfaces'};
	      for i = 1:length(possibleSurfDirNames)
		if isdir(fullfile(subjectDir,possibleSurfDirNames{i}))
		  % then return that, if it contains the WM surface file
		  val = fullfile(subjectDir,possibleSurfDirNames{i});
		  if isfile(fullfile(val,innerCoordsFilename)),return,end
		end
	      end
	      % look for the FreeSurfer directory, which should either directly contain the innerCoordsFilename
	      % or have a sub directory called surfRelax which does
	      subjectDirListing = dir(subjectDir);
	      for i = 1:length(subjectDirListing)
		if subjectDirListing(i).isdir
		  % see if the file is directly here
		  if isfile(fullfile(subjectDir,subjectDirListing(i).name,innerCoordsFilename))
		    % found it, return the directory
		    val = fullfile(subjectDir,subjectDirListing(i).name);
		    return
		  end
		  % look for surfRelaxDir
		  surfRelaxDir = fullfile(subjectDir,subjectDirListing(i).name,'surfRelax');
		  if isdir(surfRelaxDir)
		    % found surfRelaxDir
		    val = surfRelaxDir;
		    % return this one if we find the innerCoords - it is likely correct, if not
		    % will continue searching to see if we find the file in some other directory
		    if isfile(fullfile(val,innerCoordsFilename))
		      return
		    end
		  end
		end
	      end
	    end
	  end
        end
      end
    end
  case {'basecoordmap'}
    % basedata = viewGet(view,'baseCoordMap',[baseNum],[corticalDepths])
    b = getBaseNum(view,varargin);
    n = viewGet(view,'numberofbasevolumes');
    % get number of cortical depth bins
    if ieNotDefined('varargin') || (length(varargin)<2)
      corticalDepths = 0:1/(viewGet(view,'corticalDepthBins')-1):1;
    else
      corticalDepths = varargin{2};
    end
    if b & (b > 0) & (b <= n)
      val = view.baseVolumes(b).coordMap;
      % see if the coordMap is calculated for the correct number of cortical depth bins
      if ~isempty(val) && (~isfield(val,'corticalDepths') || ~isequal(val.corticalDepths,corticalDepths))
        if isfield(val,'innerCoords') && isfield(val,'outerCoords')
          % if not, then we have to do it
          %	  val.coords = (1-corticalDepth)*val.innerCoords + corticalDepth*val.outerCoords;
          val.coords = NaN([size(val.innerCoords) length(corticalDepths)]);
          cDepth=0;
          for iDepth = corticalDepths;
            cDepth=cDepth+1;
            val.coords(:,:,:,:,cDepth) = val.innerCoords + iDepth*(val.outerCoords-val.innerCoords);
          end
          val.corticalDepths = corticalDepths;
        end
      end
    end
  case {'basesurface'}
    % basedata = viewGet(view,'baseSurface',[baseNum],[corticalDepth])
    b = getBaseNum(view,varargin);
    % get cortical depth
    if ieNotDefined('varargin') || (length(varargin)<2)
      corticalDepth = viewGet(view,'corticalDepth');
    else
      corticalDepth = varargin{2};
    end
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      % get coordmap
      coordMap = view.baseVolumes(b).coordMap;
      if ~isempty(coordMap)
	if isfield(coordMap,'innerVtcs') && isfield(coordMap,'outerVtcs')
	  % find intermideate values
	  vtcs = coordMap.innerVtcs + mean(corticalDepth)*(coordMap.outerVtcs-coordMap.innerVtcs);
	  val.tris = coordMap.tris;
	elseif isfield(coordMap,'innerVtcs')
	  % if only inner is present then just return that
	  vtcs = coordMap.innerVtcs;
	  val.tris = coordMap.tris;
	else
	  vtcs = [];
	  val.tris = [];
	end
	% center surface
	if ~isempty(vtcs)
	  val.vtcs(:,1) = vtcs(:,2);
	  val.vtcs(:,2) = vtcs(:,1);
	  val.vtcs(:,3) = vtcs(:,3);
	else
	  val.vtcs = [];
	end
      else
	val = [];
      end
    end
  case {'baseclip'}
    % baseclip = viewGet(view,'baseclip',[baseNum])
    b = getBaseNum(view,varargin);
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      val = view.baseVolumes(b).clip;
    end
  case {'baserotate'}
    % baserotate = viewGet(view,'baserotate',[baseNum])
    b = getBaseNum(view,varargin);
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      baseType = viewGet(view,'baseType',b);
      baseMultiAxis = viewGet(view,'baseMultiAxis',b);
      if (baseType == 1) || ((baseType==0) && (baseMultiAxis==0))
        val = view.baseVolumes(b).rotate;
      else
	val = view.baseVolumes(b).surfaceRotate;
      end
    end
  case {'basetilt'}
    % baserotate = viewGet(view,'baseTilt',[baseNum])
    b = getBaseNum(view,varargin);
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      val = view.baseVolumes(b).tilt;
    end
  case {'basecurslice','baseslice'}
    % baseslice = viewGet(view,'baseslice',[baseNum])
    val = 1;
    b = getBaseNum(view,varargin);
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      if length(view.baseVolumes(b).curCoords) == 3
	val = view.baseVolumes(b).curCoords(viewGet(view,'baseSliceIndex',b));
      end
    end
  case {'baserawslice'}
    % baseRawSlice = viewGet(view,'baseRawSlice',orientation,sliceNum,[baseNum])
    % used to get the raw slice from the base (that is, not applying
    % any rotation that the user has set).
    if length(varargin) < 2
      disp(sprintf('(viewGet:baseRawSlice) Must specify orientation and sliceNum'));
      return
    end
    orientation = varargin{1};
    sliceNum = varargin{2};
    b = getBaseNum(view,varargin,3);
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      if any(orientation==[1 2 3])
	baseDims = size(view.baseVolumes(b).data);
	if (sliceNum>=1) && (sliceNum<=baseDims(orientation))
	  switch(orientation)
	    case {1}
	     val = squeeze(view.baseVolumes(b).data(sliceNum,:,:));
	    case {2}
	     val = squeeze(view.baseVolumes(b).data(:,sliceNum,:));
	    case {3}
	     val = squeeze(view.baseVolumes(b).data(:,:,sliceNum));
	  end
	else
	  disp(sprintf('(viewGet:baseRawSlice) sliceNum: %i out of range',sliceNum));
	end
      else
	disp(sprintf('(viewGet:baseRawSlice) orientation: %i out of range',orientation));
      end
    else
      disp(sprintf('(viewGet:baseRawSlice) baseNum: %i out of range',b));
    end
  case {'basecorticaldepth'}
    % baseslice = viewGet(view,'baseCorticalDepth',[baseNum])
    b = getBaseNum(view,varargin);
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      val = view.baseVolumes(b).curCorticalDepth;
    end
  case {'basesliceorientation'}
    % basesliceOrientation = viewGet(view,'baseSliceOrientation',[baseNum])
    b = getBaseNum(view,varargin);
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      val = view.baseVolumes(b).sliceOrientation;
    end
  case {'baserange'}
    % baserange = viewGet(view,'baserange',[baseNum])
    b = getBaseNum(view,varargin);
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      val = view.baseVolumes(b).range;
    end
  case {'basegamma'}
    % baseGamma = viewGet(view,'baseGamma',[baseNum])
    b = getBaseNum(view,varargin);
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      val = view.baseVolumes(b).gamma;
    end
  case {'baseoverlay'}
    % baseGamma = viewGet(view,'baseoverlay',[baseNum])
    b = getBaseNum(view,varargin);
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      val = view.baseVolumes(b).overlay;
    end
  case {'baseoverlayalpha'}
    % baseGamma = viewGet(view,'baseoverlayAlpha',[baseNum])
    val = 1;
    b = getBaseNum(view,varargin);
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      val = view.baseVolumes(b).overlayAlpha;
    end
  case {'basetype'}
    % baseType = viewGet(view,'baseType',[baseNum])
    % 0 is for a regular volume, 1 is for a flat, 2 is for a surface
    b = getBaseNum(view,varargin);
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      val = view.baseVolumes(b).type;
    end
  case {'basealpha'}
    % baseAlpha = viewGet(view,'baseAlpha',[baseNum])
    b = getBaseNum(view,varargin);
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      val = view.baseVolumes(b).alpha;
    end
  case {'basedisplayoverlay'}
    % baseDisplayOverlay = viewGet(view,'baseDisplayOverlay',[baseNum])
    % used for setting whether a base displays an overlay or not
    % some anatomy bases like fascicles do not always display overlays
    b = getBaseNum(view,varargin);
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      val = view.baseVolumes(b).displayOverlay;
    end
  case {'basemultidisplay'}
    % baseMultiDisplay = viewGet(view,'baseMultiDisplay',[baseNum])
    % sets whether the base will display even if it is not the
    % current basy
    b = getBaseNum(view,varargin);
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      val = view.baseVolumes(b).multiDisplay;
    end
  case {'basemultiaxis'}
    % baseMultiAxis = viewGet(view,'baseMultiDisplay',[baseNum])
    % can be 0,1 or 2 depending on whether the base should
    % display as a single, multi or 3D axis (only relevant
    % for non-flat, non-surface anatomies
    b = getBaseNum(view,varargin);
    n = viewGet(view,'numberofbasevolumes');
    if b & (b > 0) & (b <= n)
      val = view.baseVolumes(b).multiAxis;
    end
  case {'basesurfacedims'}
    % basedims = viewGet(view,'basedims',[baseNum])
    % basedims = viewGet(view,'basedims',[])
    % basedims = viewGet(view,'basedims')
    b = getBaseNum(view,varargin);
    baseType = viewGet(view,'baseType',b);
    if baseType == 2
      baseDims = view.baseVolumes(b).coordMap.dims;
      val(1) = baseDims(2);
      val(2) = baseDims(1);
      val(3) = baseDims(3);
    end
  case {'basedims'}
    % basedims = viewGet(view,'basedims',[baseNum])
    % basedims = viewGet(view,'basedims',[])
    % basedims = viewGet(view,'basedims')
    [b baseVolume] = getBaseNum(view,varargin);
    if ~isempty(baseVolume)
      val = size(baseVolume.data);
      % for an image with only one slice
      if length(val) == 2
        val(3) = 1;
      end
    end
  case {'basecurcoords'}
    % baseCurCoords = viewGet(view,'baseCurCoords',[baseNum])
    % baseCurCoords = viewGet(view,'baseCurCoords',[])
    % baseCurCoords = viewGet(view,'baseCurCoords')
    b = getBaseNum(view,varargin);
    if ~isempty(b)
      val = view.baseVolumes(b).curCoords;
    end
  case{'base2tal'}
    % xform = viewGet(view,'base2tal',[baseNum])
    % This will return the xform matrix that specifies the
    % transformation from this base to the volume
    % in talairach coordinates. Note that this also has
    % composited the shiftOriginXform.
    b = getBaseNum(view,varargin);
    n = viewGet(view,'numBase');
    if ~isempty(b) && (b > 0) && (b <= n)
      sform_code = viewGet(view,'baseSformCode',b);
      baseSform = viewGet(view,'baseSform',b);
      if ~isempty(baseSform)
        if (sform_code == 3)
          val = baseSform * shiftOriginXform;
        elseif (sform_code == 1)
          vol2mag = viewGet(view,'baseVol2mag',b);
          vol2tal = viewGet(view,'baseVol2tal',b);
          if ~isempty(vol2mag) && ~isempty(vol2tal)
            val = vol2tal * inv(vol2mag) * baseSform * shiftOriginXform;
          end
        end
      end
    end
  case{'base2mag','basexform'}
    % xform = viewGet(view,'base2mag',[baseNum])
    % The base2mag xform specifies the xform from the base
    % to the volume in magnet coordinates.
    % If the sform_code is set to 1, then this is the same
    % as baseSform *except* that the origin has been shifted
    % to start at 1,1,1. i.e. It is baseSform * shiftOriginXform
    %
    % Note that previous to adding the talairach transformation
    % code, base2mag used to be called baseXform which used
    % to be assumed to be the same as baseSform (note that
    % if you ask for basexform then shiftOriginXform is *not*
    % composited with the transformation--since this is the
    % way the old code worked).
    b = getBaseNum(view,varargin);
    n = viewGet(view,'numBase');
    if (b > 0) && (b <= n)
      sform_code = viewGet(view,'baseSformCode',b);
      baseSform = viewGet(view,'baseSform',b);
      if ~isempty(baseSform)
        if (sform_code == 1)
          val = baseSform * shiftOriginXform;
        elseif (sform_code == 3)
          vol2mag = viewGet(view,'baseVol2mag',b);
          vol2tal = viewGet(view,'baseVol2tal',b);
          if ~isempty(vol2mag) && ~isempty(vol2tal)
            val = vol2mag * inv(vol2tal) * baseSform * shiftOriginXform;
          end
        elseif (sform_code == 0)
          % If sform has not been set, then use the transform that
          % transforms this image directly on to the current anatomy
          % using the qform matrices.
          if strcmp(mrGetPref('verbose'),'Yes')
            oneTimeWarning(sprintf('noBaseSform_%i',b),...
              ['(viewGet:baseXform) sform is not set. Using qform to align '...
              'to base anatomy. Run mrAlign then mrUpdateNiftiHdr to fix this']);
          end
          baseqform = viewGet(view,'baseqform');
          val = baseqform * shiftOriginXform;
        end
      end
      % if we are being asked for baseXform
      if strcmp(lower(param),'basexform') && ~isempty(val)
        val = val * inv(shiftOriginXform);
      end
    end
  case{'base2base'}
    % xform = viewGet(view,'base2base',[baseNum1],[baseNum2])
    % This will return the xform matrix that specifies the
    % transformation from the current base coordinates to the specified
    % base (baseNum1)'s coordinates. If you specify baseNum2 it will
    % calculate the xform matrix from baseNum2 to baseNum1
    if length(varargin) < 1, disp(sprintf('(viewGet:base2base) Most specify baseNum1'));return,end
    % get the from base
    if length(varargin) == 1
      baseFromNum = viewGet(view,'curBase');
    else
      baseFromNum = varargin{2};
    end
    % make sure we have base numbers (and not names)
    baseToNum = viewGet(view,'baseNum',varargin{1});
    baseFromNum = viewGet(view,'baseNum',baseFromNum);
    % go through the magnet coordinates
    base2magTo = viewGet(view,'base2mag',baseToNum);
    base2magFrom = viewGet(view,'base2mag',baseFromNum);
    % if neither is empty, then return it
    if ~isempty(base2magTo) && ~isempty(base2magFrom)
      val = inv(base2magTo)*base2magFrom;
    else
      disp(sprintf('(viewGet:base2base) Could not compute transform of base %s to base %s',viewGet(view,'baseName',baseToNum),viewGet(view,'baseName',baseFromNum)));
    end
  case{'base2scan'}
    % xform = viewGet(view,'base2scan',[scanNum],[groupNum],[baseNum])
    % This will return the xform matrix that specifies the
    % transformation from base coordinates to scan coordinates
    % It checks whether the scan and the base are in magnet or Talairach
    % coordinates, and deals with the case when they are not in the
    % same space
    %  Note that this also has composited the shiftOriginXform.
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    b = getBaseNum(view,varargin,3);
    n = viewGet(view,'numBase');
    if (b > 0) & (b <= n) & (nscans >= s) & (s > 0)
      scan2tal = viewGet(view,'scan2tal',s,g);
      if ~isempty(scan2tal) % The scan has a Tal xform
        base2tal = viewGet(view,'base2tal',b); % check base
        if ~isempty(base2tal) % -CASE 1-: both the scan and the base have a Tal xform
          val = inv(scan2tal) * base2tal; % use it
        else %  -CASE 2-: the scan has a Tal xform but the base does not
          scan2mag = viewGet(view,'scan2mag',s,g); % check if the scan has a Mag xform
          if ~isempty(scan2mag) % if it does,
            base2mag = viewGet(view,'base2mag',b); % check if the base has a mag xform too
            if ~isempty(base2mag) % if they both do, use that, but give a warning
              oneTimeWarning(sprintf('scanBaseMismatch_%i_%i_%i',s,g,b),...
                ['WARNING: Ignoring the scan talairach transformation, '...
                'because scan '...
                'has a Talairach transformation, but the base does not. '...
                'Coordinates '...
                'are being converted through the magnet coordinate frame. '...
                'If you wanted '...
                'the transformations to use Talairach coordinates instead '...
                'of magnet '...
                'coordinates, you need to use mrAlign to export the talairach '...
                'transformation to the base']);
              val = inv(scan2mag) * base2mag;
            else % if the base doesn't have either xform, that's an error. Give a warning
              % and use the identity matrix
              oneTimeWarning(sprintf('noBaseXform_%i_%i_%i',s,g,b),...
                ['ERROR: Base does not have a transform for magnet or talairach '...
                'space. '...
                'Using the identity matrix to transform from base to scan. Run '...
                'mrAlign to fix the base.']);
              val = eye(4);
            end
          else % if the scan does not have a mag xform and the base does not have a tal xform
            oneTimeWarning(sprintf('incompatibleB1S3_%i_%i_%i',s,g,b),...
              ['Base and Scan are not compatible: Scan is in Talairach space '...
              '(and not magnet space) but Base is not. '...
              'Using the identity matrix to transform from base to scan. Run '...
              'mrAlign to get base and scan into the same space.']);
            val = eye(4);
          end
        end
      else % The scan doesn't have a Tal xform...
        scan2mag = viewGet(view,'scan2mag',s,g);
        if ~isempty(scan2mag) % ... but the scan does have a mag transform
          base2mag = viewGet(view,'base2mag',b); % check the base:
          if ~isempty(base2mag) % -CASE 3-: both base and scan have mag transform
            val = inv(scan2mag) * base2mag; % use it
            base2tal = viewGet(view,'base2tal',b); % but check if base had a Tal xform
            % so can warn user that it's being ignored
            if ~isempty(base2tal)
              oneTimeWarning(sprintf('scanBaseMismatch_%i_%i_%i',s,g,b),...
                ['WARNING: Ignoring the base talairach transformation, because '...
                'the base has a talairach '...
                'transformation but the scan does not. Coordinates are being '...
                'converted through the magnet '...
                'coordinate frame. If you want convert using the talairach '...
                'transformations, you need '...
                'to export a talairach transformation to the scan by running '...
                ' mrAlign.']);
            end
          else % -CASE 4-: Scan has a mag xform but base does not (and scan doesn't have
            % a tal xform, bc already checked that)
            oneTimeWarning(sprintf('incompatibleB3S1_%i_%i_%i',s,g,b),...
              ['Base and Scan are not compatible: Scan is in magnet space '...
              'and Base is not. Using the '...
              'qform to transform from base to scan. Run mrAlign '...
              'to get base and scan into the same space.']);
	    % Using the qforms for alignment. Note that this logic probably
	    % can be put in other sections of the base2scan where the code
	    % is returning an eye(4) matrix, but only tested here, so leaving
	    % it out in other places.
	    scanQform = viewGet(view,'scanQform');
	    baseQform = viewGet(view,'baseQform');
	    if ~isempty(scanQform) && ~isempty(baseQform)
	      val = inv(scanQform * shiftOriginXform) * baseQform * shiftOriginXform;
	    else
	      val = eye(4);
	    end
          end
        else % error if scan has neither a base nor a tal transform
          oneTimeWarning(sprintf('unknownSformCode_%i_%i_%i',s,g,b),...
            ['Scan is neither in Magnet space nor in Talairach Space.'...
            ' Using the identity matrix '...
            'to transform from base to scan. Run mrAlign to fix the scan.']);
          val = eye(4);
        end
      end
    end
  case{'base2roi'}
    % xform = viewGet(view,'base2roi',[roiNum],[baseNum])
    % This will return the xform matrix that specifies the
    % transformation from base coordinates to ROI coordinates
    % It checks whether the ROI and the base are in magnet or Talairach
    % coordinates, and deals with the case when they are not in the
    % same space
    %  Note that this also has composited the shiftOriginXform.
    b = getBaseNum(view,varargin,2);
    n = viewGet(view,'numBase');
    r = getRoiNum(view,varargin);
    nRois = viewGet(view,'numrois');
    if (b > 0) & (b <= n) & (nRois >= r) & (r > 0)
      roi2tal = viewGet(view,'roi2tal',r);
      if ~isempty(roi2tal) % The roi has a Tal xform
        base2tal = viewGet(view,'base2tal',b); % check base
        if ~isempty(base2tal) % -CASE 1-: both the roi and the base have a Tal xform
          val = inv(roi2tal) * base2tal; % use it
        else %  -CASE 2-: the roi has a Tal xform but the base does not
          roi2mag = viewGet(view,'roi2mag',r); % check if the roi has a Mag xform
          if ~isempty(roi2mag) % if it does,
            base2mag = viewGet(view,'base2mag',b); % check if the base has a mag xform too
            if ~isempty(base2mag) % if they both do, use that, but give a warning
              oneTimeWarning(sprintf('roiBaseMismatch_%i_%i',r,b),...
                ['WARNING: Ignoring the roi talairach transformation, '...
                'because roi has a Talairach transformation, but '...
                'the base does not. Coordinates are being converted '...
                'through the magnet coordinate frame. If you wanted '...
                'the transformations to use Talairach coordinates '...
                'instead of magnet coordinates, you need to use '...
                'mrAlign to export the talairach transformation to the base']);
              val = inv(roi2mag) * base2mag;
            else % if the base doesn't have either xform, that's an error.
              % Give a warning and use the identity matrix
              oneTimeWarning(sprintf('noBaseXform_%i_%i',r,b),...
                ['ERROR: Base does not have a transform for '...
                'magnet or talairach space. '...
                'Using the identity matrix to transform from '...
                'base to ROI. Run mrAlign to fix the base.']);
              val = eye(4);
            end
          else % if the roi does not have a mag xform and the base does not have a tal xform
            oneTimeWarning(sprintf('incompatibleRoiBase_%i_%i',r,b),...
              ['Base and ROI are not compatible: ROI is in Talairach space '...
              '(and not magnet space) but Base is not. '...
              'Using the identity matrix to transform from base to ROI. '...
              'Run mrAlign to get base and ROI into the same space.']);
            val = eye(4);
          end
        end
      else % The ROI doesn't have a Tal xform...
        roi2mag = viewGet(view,'roi2mag',r);
        if ~isempty(roi2mag) % ... but the ROI does have a mag transform
          base2mag = viewGet(view,'base2mag',b); % check the base:
          if ~isempty(base2mag) % -CASE 3-: both base and ROI have mag transform
            val = inv(roi2mag) * base2mag; % use it
            base2tal = viewGet(view,'base2tal',b); % but check if base had a Tal xform
            % so can warn user that it's being ignored
            if ~isempty(base2tal)
              oneTimeWarning(sprintf('roiBaseMismatch_%i_%i',r,b),...
                ['WARNING: Ignoring the base talairach transformation, '...
                'because the base has a talairach '...
                'transformation but the ROI does not. Coordinates '...
                'are being converted through the magnet '...
                'coordinate frame. If you want convert using the '...
                'talairach transformations, you need '...
                'to export a talairach transformation to the '...
                'ROI by running  mrAlign.']);
            end
          else % -CASE 4-: ROI has a mag xform but base does not (and ROI doesn't have
            % a tal xform, bc already checked that)
            oneTimeWarning(sprintf('incompatibleRoiBase_%i_%i',r,b),...
              ['Base and ROI are not compatible: ROI is in magnet '...
              'space and Base is not. Using the '...
              'identity matrix to transform from base to ROI. '...
              'Run mrAlign to get base and ROI into the same space.']);
            val = eye(4);
          end
        else % error if ROI has neither a magnet nor a tal transform
          oneTimeWarning(sprintf('unknownSformCode_%i_%i',r,b),...
            ['ROI is neither in Magnet space nor in Talairach '...
            'Space. Using the identity matrix '...
            'to transform from base to ROI. Run mrAlign to fix the ROI.']);
          val = eye(4);
        end
      end
    end
  case {'basevol2tal'}
    % basexform = viewGet(view,'baseVol2tal',[baseNum])
    % This will return the xform matrix that specifies the
    % transformation from volume coordinates to talairach
    % coordinates of the base volume that this base was aligned to.
    [b baseVolume] = getBaseNum(view,varargin);
    if ~isempty(baseVolume)
      val = baseVolume.vol2tal;
    end
  case {'basevol2mag'}
    % basexform = viewGet(view,'baseVol2mag',[baseNum])
    % This will return the xform matrix that specifies the
    % transformation from volume coordinates to magnet
    % coordinates of the base volume that this base was aligned to.
    [b baseVolume] = getBaseNum(view,varargin);
    if ~isempty(baseVolume)
      val = baseVolume.vol2mag;
    end
  case {'basesform'}
    % baseSform = viewGet(view,'baseSform',[baseNum])
    % This returns the sform in the nifit header
    % which if the sform_code is set to 1 by mrAlign
    % should be the xformation of the base to the
    % volume in magnet coordinates. Note that if
    % this base is the "canonical" volume then
    % this is the same as the qform.
    % if sform_code is set to 3 then this is the
    % trasnformation of the base to the volume in
    % talairach coordinates
    [b baseVolume] = getBaseNum(view,varargin);
    if ~isempty(baseVolume)
      val = baseVolume.hdr.sform44;
    else
      % placeholder when base images not loaded
      val = eye(4);
    end
  case {'baseqform'}
    % basexform = viewGet(view,'baseqform',[baseNum])
    [b baseVolume] = getBaseNum(view,varargin);
    if ~isempty(baseVolume)
      val = baseVolume.hdr.qform44;
    else
      % placeholder when base images not loaded
      val = eye(4);
    end
  case {'basesformcode'}
    % baseSformCode = viewGet(view,'baseSformCode',[baseNum])
    [b baseVolume] = getBaseNum(view,varargin);
    if ~isempty(baseVolume)
      val = baseVolume.hdr.sform_code;
    else
      % placeholder when base images not loaded
      val = [];
    end
  case {'basevolpermutation'}
    % basevolpermutation = viewGet(view,'basevolpermutation',[baseNum])
    [b baseVolume] = getBaseNum(view,varargin);
    if ~isempty(baseVolume)
      val = baseVolume.permutationMatrix;
    else
      % placeholder when base images not loaded
      val = eye(3);
    end
  case {'basesliceindex'}
    % basesliceindex = viewGet(view,'basesliceindex',[baseNum])
    % This admittedly arcane logic is also used in mrAlignGUI. If you
    % change this code, please make corresponding changes in that
    % function. Permutation matrix is set by loadAnat using even more
    % arcane logic.
    b = getBaseNum(view,varargin);
    baseType = viewGet(view,'baseType',b);
    % surface always has to be 1
    if baseType==2
      val = 1;
      return
    % flats have to be 3
    elseif baseType == 1
      val = 3;
      return
    end
    sliceOrientation = viewGet(view,'sliceOrientation');
    permutation = viewGet(view,'baseVolPermutation',b);
    switch sliceOrientation
      case 1   % Sagittal
        [m,index] = max(permutation * [1 0 0]');
      case 2   % Coronal
        [m,index] = max(permutation * [0 1 0]');
      case 3   % Axial
        [m,index] = max(permutation * [0 0 1]');
    end
    val = index;
  case {'basevoxelsize'}
    % basevoxelsize = viewGet(view,'basevoxelsize',[baseNum])
    [b baseVolume] = getBaseNum(view,varargin);
    if ~isempty(baseVolume)
      val = baseVolume.hdr.pixdim([2,3,4])';;
    end
    
    % ROI
  case {'visiblerois'}
    % roiList = viewGet(view,'visiblerois')
    % returns the number of all visible ROIs
    selectedROI = viewGet(view,'currentroi');
    n = viewGet(view,'numberOfROIs');
    option = viewGet(view,'showROIs');
    switch option
      case{'all','all perimeter'}
        if selectedROI
          val = [1:selectedROI-1,selectedROI+1:n,selectedROI];
        else
          val = 1:n;
        end
      case{'selected','selected perimeter'}
        val = selectedROI;
      case{'group','group perimeter'}
        val = viewGet(view,'roiGroup');
      otherwise
        val = [];
    end
    
  case {'showrois'}
    % show = viewGet(view,'showROIs')
    val = view.showROIs;
  case {'labelrois'}
    % show = viewGet(view,'showROIs')
    % returns a boolean that is set to 1 if
    % ROI labels are turned on, or 0 if not.
    val = view.labelROIs;
  case{'numberofrois','numrois','nrois'}
    % n = viewGet(view,'numberofROIs')
    val = length(view.ROIs);
  case{'currentroi','currentroinum','curroi'}
    % roiNum = viewGet(view,'currentROI')
    val = view.curROI;
  case{'roigroup','currentroigroup'}
    % roiNum = viewGet(view,'roigroup')
    % gets the roiNums of the current roi group.
    roiGroupNames = view.roiGroup;
    val = [];
    for i = 1:viewGet(view,'nrois');
      roiName = viewGet(view,'roiName',i);
      if any(strcmp(roiName,roiGroupNames))
	val(end+1) = i;
      end
    end
  case{'roigroupnames','currentroigroupnames'}
    % roiNum = viewGet(view,'roigroup')
    % gets the roiNames of the current roi group.
    if isnumeric(view.roiGroup)
      % roiGroupNames should be names not numbers (not
      % sure how this gets wrong, but fixing here
      val = {};
      for i = 1:length(view.roiGroup)
        roiName = viewGet(view,'roiName',view.roiGroup(i));
        if ~isempty(roiName)
          val{end+1} = roiName;
        end
      end
      % fix it back in the view
      view.roiGroup = val;
    else	
      val = view.roiGroup;
    end
  case{'roinum'}
    % roiNum = viewGet(view,'roiNum',roiName)
    % This will return the roiNum that correspondes
    % to the passed in name, if that roi is loaded into
    % the view. If roiName is a number, it will return
    % that number if roiNum is a valid roiNum.
    if ieNotDefined('varargin')
      mrErrorDlg('viewGet roiNum: must specify ROI name');
    end
    ROIname = varargin{1};
    % if no ROIs are loaded than return empty
    nROIs = viewGet(view,'nROIs');
    if nROIs == 0
      return
    end
    % check for number
    if isscalar(ROIname)
      if (ROIname >= 1) && (ROIname <= viewGet(view,'nROIs'))
        val = round(ROIname);
      end
      return
    end
    ROInames = {view.ROIs(:).name};
    % note that it is possible to have more than one
    % ROI with the same name, so pick the last one in
    % the list, since we usually are getting the number
    % for the last ROI that was just created (This should not
    % happen anymore since viewSet checks for duplicates).
    val = find(strcmp(ROIname,ROInames));
    if length(val) > 1
      mrWarnDlg(sprintf('(viewGet) you have multiple ROIs with the name %s',ROIname));
      val = last(val);
    end
  case{'roi'}
    % roi = viewGet(view,'roi',[roiNum])
    r = getRoiNum(view,varargin);
    n = viewGet(view,'numberofROIs');
    if r & (r > 0) & (r <= n)
      val = view.ROIs(r);
    end
  case {'roinames'}
    % roiNames = viewGet(view,'roiNames')
    if ieNotDefined('varargin')
      roiNum = viewGet(view,'currentROI');
    else
      roiNum = varargin{1};
    end
    if isempty(roiNum)
      roiNum = viewGet(view,'currentROI');
    end
    if ~isempty(view.ROIs)
      val = {view.ROIs(:).name};
    end
  case{'roiname'}
    % roiName = viewGet(view,'roiName',[roiNum])
    r = getRoiNum(view,varargin);
    % handle cell arrays
    if iscell(r)
      for i = 1:length(r)
        val(end+1) = viewGet(view,'roiName',r{i});
      end
      return
    end
    % if a string, just validate that it really is an roi
    if isstr(r)
      if ~isempty(viewGet(view,'roiNum',r))
        val = r;
      else
        val = [];
      end
      return
    end
    n = viewGet(view,'numberofROIs');
    if r & (r > 0) & (r <= n)
      if length(r) == 1
	val = view.ROIs(r).name;
      else
	val = {view.ROIs(r).name};
      end
    end
  case{'roicoords'}
    % roiCoords = viewGet(view,'roiCoords',[roiNum])
    r = getRoiNum(view,varargin);
    n = viewGet(view,'numberofROIs');
    if r & (r > 0) & (r <= n)
      val = view.ROIs(r).coords;
    end
  case{'prevroicoords'}
    % prevRoiCoords = viewGet(view,'prevRoiCoords')
    val = view.prevROIcoords;
  case{'roicolor'}
    % roicolor = viewGet(view,'roicolor',[roiNum])
    r = getRoiNum(view,varargin);
    n = viewGet(view,'numberofROIs');
    if r & (r > 0) & (r <= n)
      val = view.ROIs(r).color;
    end
  case{'roicolorrgb'}
    % roicolor = viewGet(view,'roicolorRGB',[roiNum])
    % roicolor = viewGet(view,'roicolorRGB',[])
    % roicolor = viewGet(view,'roicolorRGB')
    if ieNotDefined('varargin')
      roicolor = viewGet(view,'roiColor');
    else
      roicolor = viewGet(view,'roiColor',varargin{1});
    end
    val = color2RGB(roicolor);
  case{'roivol2mag'}
    % roiVol2mag = viewGet(view,'roiVol2mag',[roiNum])
    % returns the xform of the volume coordinates to
    % the magnet coordinates for the canonical volume
    r = getRoiNum(view,varargin);
    n = viewGet(view,'numberofROIs');
    if r & (r > 0) & (r <= n)
      val = view.ROIs(r).vol2mag;
    end
  case{'roivol2tal'}
    % roiVol2tal = viewGet(view,'roiVol2tal',[roiNum])
    % returns the xform of the volume coordinates to
    % the talairach coordinates for the canonical volume
    r = getRoiNum(view,varargin);
    n = viewGet(view,'numberofROIs');
    if r & (r > 0) & (r <= n)
      val = view.ROIs(r).vol2tal;
    end
  case{'roisformcode'}
    % roiSformCode = viewGet(view,'roiSformCode',[roiNum])
    % returns the sFormCode of the transform saved in
    % view.ROIs(r).sformcode, which is the sformcode of the
    % base on which the ROI was originally defined.
    r = getRoiNum(view,varargin);
    n = viewGet(view,'numberofROIs');
    if r & (r > 0) & (r <= n)
      val = view.ROIs(r).sformCode;
    end
  case{'roi2mag','roixform'}
    % roi2mag = viewGet(view,'roi2mag',[roiNum])
    % roi2mag returns the transformation of the roi
    % to the volume in magnet coordinates. It has
    % composited on it the shiftOriginXform.
    % The only difference between this and roiXform
    % is the incluson of the shiftOriginXform
    r = getRoiNum(view,varargin);
    n = viewGet(view,'numberofROIs');
    if r & (r > 0) & (r <= n)
      sform_code = viewGet(view,'roiSformCode',r);
      if (sform_code == 1)
        val = view.ROIs(r).xform * shiftOriginXform;
        % if the sform_code was set to 0, then either
        % it has a xform that came from the bases qform
        % or none at all (in which case return identity)
      elseif (sform_code == 0)
        if isempty(view.ROIs(r).xform)
          val = eye(4) * shiftOriginXform;
        else
          val = view.ROIs(r).xform * shiftOriginXform;
        end
      elseif (sform_code == 3)
        vol2tal = viewGet(view,'roiVol2tal',r);
        vol2mag = viewGet(view,'roiVol2mag',r);
        if ~isempty(vol2tal) && ~isempty(vol2mag)
          val = vol2mag * inv(vol2tal) * view.ROIs(r).xform * shiftOriginXform;
        end
      end
      if strcmp(lower(param),'roixform') && ~isempty(val)
        val = val * inv(shiftOriginXform);
      end
    else
      val = eye(4);
    end
  case{'roi2tal'}
    % roi2tal = viewGet(view,'roi2tal',[roiNum])
    % roi2tal returns the transformation of the roi
    % to the volume in talairach coordinates. It has
    % composited on it the shiftOriginXform.
    r = getRoiNum(view,varargin);
    n = viewGet(view,'numberofROIs');
    if r & (r > 0) & (r <= n)
      sform_code = viewGet(view,'roisformCode',r);
      if (sform_code == 3)
        val = view.ROIs(r).xform * shiftOriginXform;
      elseif (sform_code == 1)
        vol2tal = viewGet(view,'roiVol2tal',r);
        vol2mag = viewGet(view,'roiVol2mag',r);
        if ~isempty(vol2tal) && ~isempty(vol2mag)
          val = vol2tal * inv(vol2mag) * view.ROIs(r).xform * shiftOriginXform;
        end
      end
    else
      val = eye(4);
    end
  case{'roi2roi'}
    % xform = viewGet(view,'roi2roi',[fromROINum],[toROINum])
    % This will return the xform matrix that specifies the
    % transformation from coordinates of the 'From' ROI to coordinates of the 'To' ROI
    % (E.g., given x,y,s from the 'From' ROI, multiply by the xform calculated in this
    % call to convert to x',y',s' in the 'To' ROI.)
    % It checks whether the ROIs are each in magnet or Talairach
    % coordinates, and deals with the case when they are not in the
    % same space
    %  Note that this also has composited the shiftOriginXform.
    rFrom = getRoiNum(view,varargin);
    rTo = getRoiNum(view,varargin,2);
    n = viewGet(view,'numberofROIs');
    if (n >= rFrom) && (rFrom > 0) && (n >= rTo) && (rTo > 0)
      roi2talFrom = viewGet(view,'roi2tal',rFrom); % check if FROMroi has talXform
      if ~isempty(roi2talFrom) % if the FROMroi has a Tal xform
        roi2talTo = viewGet(view,'roi2tal',rTo); % check TOroi
        if ~isempty(roi2talTo) % -CASE 1-: both rois have a Tal xform, then use it
          % first check if they are the same, then just
          % return identity matrix, so there's no roundoff
          % error when compute the inverse
          if isequal(roi2talTo, roi2talFrom)
            val = eye(4);
          else
            val = inv(roi2talTo) * roi2talFrom;
          end
        else %  -CASE 2-: the FROM roi has a Tal xform but the TO roi doesn't
          roi2magFrom = viewGet(view,'roi2mag',rFrom); % check if the FROMroi has magXform
          if ~isempty(roi2magFrom) % if it does,
            roi2magTo = viewGet(view,'roi2mag',rTo); % check if the TOroi has magXform
            if ~isempty(roi2magTo) % if they both do, use that, but give a warning
              oneTimeWarning(sprintf('roiMismatch_%i_%i',rFrom,rTo),...
                ['WARNING: Ignoring the ROI talairach transformation, '...
                'because FROM ROI has a Talairach transformation,'...
                ' but TO ROI does not. Coordinates are being converted'...
                ' through the magnet coordinate frame.'...
                ' If you wanted the transformations to use Talairach'...
                ' coordinates instead of magnet coordinates, you'...
                ' need to use mrAlign to export the talairach '...
                ' transformation to the TO ROI']);
              if isequal(roi2magTo,roi2magFrom) % if they're the same, avoid
                % the roundoff errors of calling inv
                val = eye(4);
              else
                val = inv(roi2magTo) * roi2magFrom;
              end
              
            else % if the TOroi doesn't have either xform, that's an error.
              % Give a warning and use the identity matrix
              oneTimeWarning(sprintf('noRoiXform_%i_%i',rFrom,rTo),...
                ['ERROR: TO ROI does not have a transform for magnet '...
                'or talairach space. Using the identity '...
                'matrix to transform from ROI to ROI. '...
                'Run mrAlign to fix the TO ROI.']);
              val = eye(4);
            end
          else % if the FROM ROI does not have a mag xform
            % and the TO ROI does not have a tal xform
            oneTimeWarning(sprintf('incompatibleROIs_%i_%i',rFrom,rTo),...
              ['ROIs are not compatible: FROM ROI is in '...
              'Talairach space but TO ROI is not. '...
              'Using the identity matrix to transform from ROI to ROI. '...
              'Run mrAlign to get ROIs into the same space.']);
            val = eye(4);
          end
        end
      else % The FROM ROI doesn't have a Tal xform...
        roi2magFrom = viewGet(view,'roi2mag',rFrom);
        if ~isempty(roi2magFrom) % ... but the FROM ROI does have a mag transform
          roi2magTo = viewGet(view,'roi2mag',rTo); % check the TO ROI
          if ~isempty(roi2magTo) % -CASE 3-: both ROIs have mag transform
            if isequal(roi2magTo,roi2magFrom) % if they're the same, avoid
              % the roundoff errors of calling inv
              val = eye(4);
            else
              val = inv(roi2magTo) * roi2magFrom;
            end
            % but check if TO scan had a Tal xform so can warn user that it's being ignored
            roi2talTo = viewGet(view,'roi2tal',rTo);
            if ~isempty(roi2talTo)
              oneTimeWarning(sprintf('roiMismatch_%i_%i',rFrom,rTo),...
                ['WARNING: Ignoring the TO ROIs talairach transformation, '...
                'because the TO ROI has a talairach transformation but '...
                'the FROM ROI does not. Coordinates are being converted '...
                'through the magnet coordinate frame. '...
                'If you want convert using the talairach transformations, '...
                'you need to export a talairach transformation to '...
                'the TO ROI by running  mrAlign.']);
            end
          else % -CASE 4-: FROM ROI has a mag xform but TO ROI does not
            % (and FROM ROI doesn't have a tal xform, bc already checked that)
            oneTimeWarning(sprintf('incompatibleROIs_%i_%i',rFrom,rTo),...
              ['ROIs are not compatible: FROM ROI is in magnet space but '...
              'TO ROI is not. Using the identity matrix to transform from '...
              'ROI to ROI. Run mrAlign to get ROIs into the same space.']);
            val = eye(4);
          end
        else % error if ROI has neither a base nor a tal transform
          oneTimeWarning(sprintf('unknownSformCode_%i_%i',rFrom,rTo),...
            ['FROM ROI is neither in Magnet space nor in Talairach Space.'...
            ' Using the identity matrix to transform from ROI to ROI.'...
            ' Run mrAlign to fix the ROI.']);
          val = eye(4);
        end
      end
    end
  case{'roinotes'}
    % roiNotes = viewGet(view,'roinotes',[roiNum])
    r = getRoiNum(view,varargin);
    n = viewGet(view,'numberofROIs');
    if r & (r > 0) & (r <= n)
      val = view.ROIs(r).notes;
    else
      val = '';
    end
  case{'roicreatedby'}
    % roiCreatedBy = viewGet(view,'roiCreatedBy',[roiNum])
    r = getRoiNum(view,varargin);
    n = viewGet(view,'numberofROIs');
    if r & (r > 0) & (r <= n)
      val = view.ROIs(r).createdBy;
    else
      val = '';
    end
  case{'roicreatedonbase'}
    % roiCreatedOnBase = viewGet(view,'roiCreatedOnBase',[roiNum])
    r = getRoiNum(view,varargin);
    n = viewGet(view,'numberofROIs');
    if r & (r > 0) & (r <= n)
      val = view.ROIs(r).createdOnBase;
    else
      val = '';
    end
  case{'roidisplayonbase'}
    % roiDisplayOnBase = viewGet(view,'roiDisplayOnBase',[roiNum])
    r = getRoiNum(view,varargin);
    n = viewGet(view,'numberofROIs');
    if r & (r > 0) & (r <= n)
      val = view.ROIs(r).displayOnBase;
    else
      val = '';
    end
  case{'roicreatedfromsession'}
    % roiCreatedFromSession = viewGet(view,'roiCreatedFromSession',[roiNum])
    r = getRoiNum(view,varargin);
    n = viewGet(view,'numberofROIs');
    if r & (r > 0) & (r <= n)
      val = view.ROIs(r).createdFromSession;
    else
      val = '';
    end
  case{'roibranchnum'}
    % roiBranchNum = viewGet(view,'roiBranchNum',[roiNum])
    r = getRoiNum(view,varargin);
    n = viewGet(view,'numberofROIs');
    if r & (r > 0) & (r <= n)
      val = view.ROIs(r).branchNum;
    else
      val = '';
    end
  case{'roisubjectid'}
    % roiSubjectID = viewGet(view,'roiSubjectID',[roiNum])
    r = getRoiNum(view,varargin);
    n = viewGet(view,'numberofROIs');
    if r & (r > 0) & (r <= n)
      val = view.ROIs(r).subjectID;
    else
      val = '';
    end
  case{'roivoxelsize'}
    % roivoxelsize = viewGet(view,'roivoxelsize',[roiNum])
    r = getRoiNum(view,varargin);
    n = viewGet(view,'numberofROIs');
    if r & (r > 0) & (r <= n)
      val = view.ROIs(r).voxelSize;
    end
  case('roivolume')
    % roivolume = viewGet(view,'roivolume',[roiNum])
    r = getRoiNum(view,varargin);
    n = viewGet(view,'numberofROIs');
    if r & (r > 0) & (r <= n)
      roiVoxelSize = view.ROIs(r).voxelSize;
      roiCoords = view.ROIs(r).coords;
      if ~isempty(roiVoxelSize) & ~isempty(roiCoords)
        val = prod(roiVoxelSize) * size(roiCoords,2);
      end
    end
  case{'roidate'}
    % roidate = viewGet(view,'roidate',[roiNum])
    r = getRoiNum(view,varargin);
    n = viewGet(view,'numberofROIs');
    if r & (r > 0) & (r <= n)
      val = view.ROIs(r).date;
    end
    
    % Cache
  case{'roicacheid'}
    % cacheID = viewGet(view,'ROICacheID')
    % cacheID = viewGet(view,'ROICacheID',roinum)
    % cacheID = viewGet(view,'ROICacheID',roinum,basenum)
    if length(varargin)<1
      % if we are not passed in a specific ROI, then
      % we are doing all ROIs
      roiNums = 1:length(view.ROIs);
    else
      roiNums = varargin{1};
    end
    if length(varargin)<2
      baseNum = viewGet(view,'curBase');
    else
      baseNum = varargin{2};
    end
    rotate = viewGet(view,'rotate');
    % go through all the bases and look for
    % the first one with a matching xform and
    % voxel size. We do this because flat maps
    % all have the same baseXform and voxle size
    % and so there is no need to recompute them
    % for each base map
    baseXform = viewGet(view,'baseXform',baseNum);
    baseVoxelSize = viewGet(view,'baseVoxelSize',baseNum);
    for bnum = 1:viewGet(view,'numBase')
      if (isequal(baseXform,viewGet(view,'baseXform',bnum)) && ...
          isequal(baseVoxelSize,viewGet(view,'baseVoxelSize',bnum)))
        baseMatch = bnum;
        break;
      end
    end
    baseName = viewGet(view,'baseName',baseMatch);
    val = sprintf('%s_%i',baseName,rotate);
    for i = roiNums
      val = sprintf('%s_%s_%i',val,view.ROIs(i).name,size(view.ROIs(i).coords,2));
    end
  case{'roicache'}
    % cacheVal = viewGet(view,'ROICache')
    % cacheVal = viewGet(view,'ROICache',roinum)
    % cacheVal = viewGet(view,'ROICache',roinum,basenum)
    roiID = viewGet(view,'ROICacheID',varargin{:});
    % and retrieve from the cache
    [val MLR.caches{view.viewNum}.roiCache] = ...
      mrCache('find',MLR.caches{view.viewNum}.roiCache,roiID);
  case{'basecacheid'}
    % cacheID = viewGet(view,'baseCacheID')
    % cacheID = viewGet(view,'baseCacheID',baseNum);
    % cacheID = viewGet(view,'baseCacheID',baseNum,sliceNum);
    % cacheID = viewGet(view,'baseCacheID',baseNum,sliceNum,sliceIndex);
    if length(varargin)<1
      baseNum = viewGet(view,'curBase');
    else
      baseNum = varargin{1};
    end
    if length(varargin)>=2
      currentSlice = varargin{2};
    else
      currentSlice = viewGet(view,'curSlice');
    end
    if length(varargin)>=3
      sliceIndex = varargin{3};
    else
      sliceIndex = viewGet(view,'baseSliceIndex');
    end
    if length(varargin)>=4
      rotate = varargin{4};
    else
      rotate = viewGet(view,'rotate');
    end
    baseName = viewGet(view,'baseName',baseNum);
    clip = viewGet(view,'baseClip',baseNum);
    gamma = viewGet(view,'baseGamma',baseNum);
    % only use the corticalDepth if this is a flat
    if viewGet(view,'baseType')
      corticalDepthBins = viewGet(view,'corticalDepthBins');
    else
      corticalDepthBins = 0;
    end
    val = sprintf('%s_%i_%i_%i_%s_%s_%0.2f_%0.2f',baseName,currentSlice,sliceIndex,rotate,num2str(clip(1)),num2str(clip(end)),gamma,corticalDepthBins);
  case{'basecache'}
    % cacheVal = viewGet(view,'baseCache')
    % cacheVal = viewGet(view,'baseCache',baseNum)
    % cacheVal = viewGet(view,'baseCache',baseNum,sliceNum)
    % cacheVal = viewGet(view,'baseCache',baseNum,sliceNum,sliceIndex)
    baseID = viewGet(view,'baseCacheID',varargin{:});
    % and retrieve from the cache
    [val MLR.caches{view.viewNum}.baseCache] = ...
      mrCache('find',MLR.caches{view.viewNum}.baseCache,baseID);
  case{'overlaycacheid'}
    % cacheID = viewGet(view,'overlayCacheID')
    % cacheID = viewGet(view,'overlayCacheID',baseNum)
    % cacheID = viewGet(view,'overlayCacheID',baseNum,sliceNum)
    % cacheID = viewGet(view,'overlayCacheID',baseNum,sliceNum,sliceIndex)
    % cacheID = viewGet(view,'overlayCacheID',baseNum,sliceNum,sliceIndex,rotate)
    %curSlice = viewGet(view,'curSlice');
    %analysisNum = viewGet(view,'currentAnalysis');
    %curOverlay = viewGet(view,'currentOverlay');
    %clip = viewGet(view,'overlayClip',curOverlay);
    %overlayType = viewGet(view,'overlayCtype',curOverlay);
    %overlayRange = viewGet(view,'overlayRange',curOverlay);
    % forgoe viewgets, and just grab stuff here explicitly
    % this saves about 100ms
    val = -1;
    analysisNum = view.curAnalysis;
    if ~isempty(analysisNum)
      if length(varargin)<1
	baseNum = viewGet(view,'curBase');
      else
	baseNum = varargin{1};
      end
      if length(varargin)>=2
	curSlice = varargin{2};
      else
	curSlice = viewGet(view,'curSlice');
      end
      if length(varargin)>=3
	sliceIndex = varargin{3};
      else
	sliceIndex = viewGet(view,'baseSliceIndex');
      end
      if length(varargin)>=4
	rotate = varargin{4};
      else
	rotate = viewGet(view,'rotate');
      end
      baseName = viewGet(view,'baseName',baseNum);
      curOverlay = view.analyses{analysisNum}.curOverlay;
      if ~isempty(curOverlay)
        % get all clips
        clip = [];
        for i = 1:length(view.analyses{analysisNum}.overlays)
          clip = [clip view.analyses{analysisNum}.overlays(i).clip];
        end
        overlayRange = view.analyses{analysisNum}.overlays(curOverlay).range;
        scanNum = viewGet(view,'curScan');
        alpha = viewGet(view,'alpha');
        alphaOverlay = char(viewGet(view,'alphaOverlay'))';
        alphaOverlayExponent = viewGet(view,'alphaOverlayExponent');
        clipAcrossOverlays = viewGet(view,'clipAcrossOverlays');
        multiSliceProjection = mrGetPref('multiSliceProjectionMethod');
        % need to recalculate overlay if this is aflat
        % and the cortical depth has changed
        if viewGet(view,'baseType',baseNum)
          corticalDepth = viewGet(view,'corticalDepth');
        else
          corticalDepth = 0;
        end
	baseOverlayAlpha = viewGet(view,'baseOverlayAlpha',baseNum);
	if length(baseOverlayAlpha)>1
	  % use the sum of all alpha as a proxy for computing cache
	  % this might fail sometimes to produce unique values, but
	  % should be ok for most situations
	  baseOverlayAlpha = sum(baseOverlayAlpha(:));
	end
	baseOverlay = viewGet(view,'baseOverlay',baseNum);
	if isempty(baseOverlay)
	  baseOverlay = '';
	elseif ~isstr(baseOverlay)
	  % if baseOverlay is actually RGB values, then take
	  % the sum - again, this could fail be to unique, but
	  % likely will be ok in most situations
	  baseOverlay = sprinf('%f',sum(baseOverlay(:)));
	end
        % calculate string
        val = sprintf('%i_%s_%i_%i_%i_%s_%s_%s_%i_%i_%s_%i_%s_%f_%s_%s_%i',scanNum,baseName,curSlice,sliceIndex,analysisNum,mat2str(curOverlay),mat2str(clip),mat2str(overlayRange),rotate,alpha,mat2str(corticalDepth),clipAcrossOverlays,multiSliceProjection,baseOverlayAlpha,baseOverlay,alphaOverlay(:)',mat2str(alphaOverlayExponent));
      end
    end
    %    val = curSlice*analysisNum*curOverlay;
  case{'overlaycache'}
    % cacheVal = viewGet(view,'overlayCache')
    % cacheVal = viewGet(view,'overlayCache',baseNum)
    % cacheVal = viewGet(view,'overlayCache',baseNum,sliceNum)
    % cacheVal = viewGet(view,'overlayCache',baseNum,sliceNum,sliceIndex)
    % cacheVal = viewGet(view,'overlayCache',baseNum,sliceNum,sliceIndex,rotate)

    % get the overlay ID
    overlayID = viewGet(view,'overlayCacheID',varargin{:});
    % and retrieve from the cache
    [val MLR.caches{view.viewNum}.overlayCache] = ...
      mrCache('find',MLR.caches{view.viewNum}.overlayCache,overlayID);
    
    
    % analysis
  case{'numberofanalyses','nanalyses','numanalyses'}
    % n = viewGet(view,'numberofAnalyses')
    % n = viewGet(view,'numberofAnalyses',[groupNum])
    if ieNotDefined('varargin')
      val = length(view.analyses);
    else
      % if this is the current group
      % then just get it
      if varargin{1} == viewGet(view,'curGroup')
        val = length(view.analyses);
      elseif (varargin{1} > 0) && (varargin{1} <= viewGet(view,'nGroups'))
        val = length(view.loadedAnalyses{varargin{1}});
      else
        val = [];
      end
    end
  case{'currentanalysis','curanalysis'}
    % n = viewGet(view,'currentAnalysis')
    val = view.curAnalysis;
  case{'analysisnum'}
    % n = viewGet(view,'analysisNum',analysisName)
    if ieNotDefined('varargin')
      mrErrorDlg('viewGet analysisNum: must specify analysisName');
    end
    analysisName = varargin{1};
    if ~isempty(view.analyses)
      analysisNames = viewGet(view,'analysisnames');
      val = find(strcmp(analysisName,analysisNames));
    end
  case {'analysis'}
    % analysis = viewGet(view,'analysis',[analysisNum])
    % analysis = viewGet(view,'analysis',[])
    % analysis = viewGet(view,'analysis')
    % analysis = viewGet(view,'analysis',[analysisNum],[groupNum])
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{1};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    % if the user passed in group num, then that
    % means to check the "loadedAnalyses" for that group
    if length(varargin) <= 1
      groupNum = viewGet(view,'curGroup');
    else
      groupNum = varargin{2};
    end
    n = viewGet(view,'numberofAnalyses',groupNum);
    if analysisNum & (analysisNum > 0) & (analysisNum <= n)
      if groupNum == viewGet(view,'curGroup')
        val = view.analyses{analysisNum};
      else
        val = view.loadedAnalyses{groupNum}{analysisNum};
      end
    end
  case {'d'}
    % The d structure is a data structure that contains
    % outcomes of analyses that are not overlays. For
    % example, event-related analyses have a d strucutre
    % with the estimated hemodynamic responses and other
    % computed fields
    % d = viewGet(view,'d',[scanNum],[analysisNum])
    if ieNotDefined('varargin')
      scanNum = viewGet(view,'currentScan');
    else
      scanNum = varargin{1};
    end
    if isempty(scanNum)
      scanNum = viewGet(view,'currentScan');
    end
    % get analysis num
    if length(varargin) <= 1
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{2};
    end
    n = viewGet(view,'numberofAnalyses');
    if analysisNum & (analysisNum > 0) & (analysisNum <= n)
      if isfield(view.analyses{analysisNum},'d')
        if ~isempty(scanNum) && (scanNum > 0 ) && (scanNum <= length(view.analyses{analysisNum}.d))
          val = view.analyses{analysisNum}.d{scanNum};
        end
      end
    end
  case {'analysisnames'}
    % analysisNames = viewGet(view,'analysisNames',[groupnum])
    if ieNotDefined('varargin')
      analyses = view.analyses;
      groupnum = viewGet(view,'currentGroup');
    elseif length(varargin)==1
      groupnum = varargin{1};
      if groupnum == viewGet(view,'currentGroup')
        analyses = view.analyses;
      else
        analyses = viewGet(view,'loadedanalyses',varargin{1});
      end
    end
    if ~isempty(analyses)
      numAnalyses = viewGet(view,'numberofAnalyses',groupnum);
      names = cell(1,numAnalyses);
      for a = 1:numAnalyses
        names{a} = analyses{a}.name;
      end
      val = names;
    end
  case {'analysisname'}
    % analysisName = viewGet(view,'analysisName',[analysisNum])
    % analysisName = viewGet(view,'analysisName',[])
    % analysisName = viewGet(view,'analysisName')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{1};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    n = viewGet(view,'numberofAnalyses');
    if analysisNum & (analysisNum > 0) & (analysisNum <= n)
      val = view.analyses{analysisNum}.name;
    end
  case {'analysistypes'}
    % analysisTypes = viewGet(view,'analysisTypes')
    if ~isempty(view.analyses)
      numAnalyses = viewGet(view,'numberofAnalyses');
      types = cell(1,numAnalyses);
      for a = 1:numAnalyses
        types{a} = view.analyses{a}.type;
      end
      val = types;
    end
  case {'analysistype'}
    % analysisType = viewGet(view,'analysisType',[analysisNum])
    % analysisType = viewGet(view,'analysisType',[])
    % analysisType = viewGet(view,'analysisType')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{1};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    n = viewGet(view,'numberofAnalyses');
    if analysisNum & (analysisNum > 0) & (analysisNum <= n)
      val = view.analyses{analysisNum}.type;
    end
  case {'analysisgroupname'}
    % analysisGroup = viewGet(view,'analysisGroupName',[analysisNum])
    % analysisGroup = viewGet(view,'analysisGroupName',[])
    % analysisGroup = viewGet(view,'analysisGroupName')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{1};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    n = viewGet(view,'numberofAnalyses');
    if analysisNum & (analysisNum > 0) & (analysisNum <= n)
      val = view.analyses{analysisNum}.groupName;
    end
  case {'analysisfunction'}
    % analysisFunction = viewGet(view,'analysisFunction',[analysisNum]);
    % analysisFunction = viewGet(view,'analysisFunction',[]);
    % analysisFunction = viewGet(view,'analysisFunction');
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{1};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    n = viewGet(view,'numberofAnalyses');
    if analysisNum & (analysisNum > 0) & (analysisNum <= n)
      val = view.analyses{analysisNum}.function;
    end
  case {'analysisreconcilefunction','reconcilefunction'}
    % reconcileFunction = viewGet(view,'reconcilefunction',[analysisNum])
    % reconcileFunction = viewGet(view,'reconcilefunction',[])
    % reconcileFunction = viewGet(view,'reconcilefunction')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{1};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    n = viewGet(view,'numberofAnalyses');
    if analysisNum & (analysisNum > 0) & (analysisNum <= n)
      val = view.analyses{analysisNum}.reconcileFunction;
    end
  case {'analysismergefunction','mergefunction'}
    % mergeFunction = viewGet(view,'mergefunction',[analysisNum])
    % mergeFunction = viewGet(view,'mergefunction',[])
    % mergeFunction = viewGet(view,'mergefunction')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{1};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    n = viewGet(view,'numberofAnalyses');
    if analysisNum & (analysisNum > 0) & (analysisNum <= n)
      val = view.analyses{analysisNum}.mergeFunction;
    end
  case {'analysisguifunction'}
    % analysisGui = viewGet(view,'analysisGuiFunction',[analysisNum])
    % analysisGui = viewGet(view,'analysisGuiFunction',[])
    % analysisGui = viewGet(view,'analysisGuiFunction')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{1};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    n = viewGet(view,'numberofAnalyses');
    if analysisNum & (analysisNum > 0) & (analysisNum <= n)
      val = view.analyses{analysisNum}.guiFunction;
    end
  case {'analysisparams'}
    % analysisParams = viewGet(view,'analysisParams',[analysisNum])
    % analysisParams = viewGet(view,'analysisParams',[])
    % analysisParams = viewGet(view,'analysisParams')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{1};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    n = viewGet(view,'numberofAnalyses');
    if analysisNum & (analysisNum > 0) & (analysisNum <= n)
      val = view.analyses{analysisNum}.params;
    end
  case {'analysisdata'}
    % analysisData = viewGet(view,'analysisData',[analysisNum])
    % analysisData = viewGet(view,'analysisData',[])
    % analysisData = viewGet(view,'analysisData')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{1};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    n = viewGet(view,'numberofAnalyses');
    if analysisNum & (analysisNum > 0) & (analysisNum <= n)
      val = view.analyses{analysisNum}.d;
    end
  case {'analysisdate'}
    % analysisdate = viewGet(view,'analysisDate',[analysisNum])
    % analysisdate = viewGet(view,'analysisDate',[])
    % analysisdate = viewGet(view,'analysisDate')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{1};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    n = viewGet(view,'numberofAnalyses');
    if analysisNum & (analysisNum > 0) & (analysisNum <= n)
      val = view.analyses{analysisNum}.date;
    end
  case {'overlayinterpfunction'}
    % fnctn = viewGet(view,'overlayInterpFunction',[analysisNum])
    % fnctn = viewGet(view,'overlayInterpFunction',[])
    % fnctn = viewGet(view,'overlayInterpFunction')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{1};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    n = viewGet(view,'numberofAnalyses');
    if analysisNum & (analysisNum > 0) & (analysisNum <= n)
      if isfield(view.analyses{analysisNum},'overlayInterpFunction');
        val = view.analyses{analysisNum}.overlayInterpFunction;
      else
        val = [];
      end
    end
    
    % overlay
  case{'overlays'}
    % overlays = viewGet(view,'overlays',[analysisNum])
    % overlays = viewGet(view,'overlays',[])
    % overlays = viewGet(view,'overlays')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{1};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis)
      val = analysis.overlays;
    end
  case{'numberofoverlays','numoverlays','noverlays'}
    % n = viewGet(view,'numberofOverlays',[analysisNum])
    % n = viewGet(view,'numberofOverlays',[])
    % n = viewGet(view,'numberofOverlays')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{1};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis)
      val = length(analysis.overlays);
    end
  case{'currentoverlay','curoverlay'}
    % n = viewGet(view,'currentOverlay',[analysisNum])
    % n = viewGet(view,'currentOverlay',[])
    % n = viewGet(view,'currentOverlay')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{1};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis)
      val = analysis.curOverlay;
    end
    
  case {'curclippingoverlay','currentclippingoverlay'}
    % curclippingoverlay = viewGet(view,'curclippingoverlay')
    if ~isempty(viewGet(view,'fignum'))
      handles = guidata(view.figure);
    else
      handles=[];
    end
    if isempty(handles) || ~isfield(handles,'clippingOverlaysListbox') 
      val=viewGet(view,'currentOverlay');
    else
      clippingOverlays = get(handles.clippingOverlaysListbox,'string');
      % check for no current overlays
      if isempty(clippingOverlays)
	val = [];
      else
	clippingOverlayNum = get(handles.clippingOverlaysListbox,'value');
	if isempty(clippingOverlayNum) %this can happen on a Mac
	  clippingOverlayNum = 1;
	  set(handles.clippingOverlaysListbox,'value',1);
	end
	clippingOverlay = clippingOverlays{clippingOverlayNum};
	val = viewGet(view,'overlayNum',clippingOverlay(3:end)); %ignore the first 2 characters which are used to put a star before the name (see mlrGuiSet, case 'clippingOverlays'))
      end
    end
    
  case {'clippingoverlaylist'}
    if ~isempty(viewGet(view,'fignum'))
      handles = guidata(view.figure);
    else
      handles=[];
    end
    if ~isempty(handles) && isfield(handles,'clipAcrossOverlays') && get(handles.clipAcrossOverlays,'value')
      clippingOverlayList = 1:viewGet(view,'nOverlays');
    else
      %set overlay and alpha overlay names in clipping box 
      curOverlays = viewGet(view,'curOverlay');
      analysisNum = viewGet(view,'currentAnalysis');
      clippingOverlayList = curOverlays;
      for iOverlay=1:length(curOverlays)
        alphaOverlayNum = viewGet(view,'overlayNum',viewGet(view,'alphaOverlay',curOverlays(iOverlay),analysisNum),analysisNum);
        if ~isempty(alphaOverlayNum)
          clippingOverlayList(end+1)=alphaOverlayNum;
        end
      end
    end
    %unique without ordering
    [dummy,uniqueIndices]=unique(clippingOverlayList,'first');
    val=clippingOverlayList(unique(uniqueIndices));
  
  case{'overlaynum'}
    % n = viewGet(view,'overlayNum',overlayName(s),[analysisNum])
    % n = viewGet(view,'overlayNum',overlayName(s),[])
    % n = viewGet(view,'overlayNum',overlayName(s))
    %     overlayName(s) may be a string or cell array of strings
    if ieNotDefined('varargin')
      mrErrorDlg('viewGet overlayNum: must specify overlayName.');
    end
    overlayName = varargin{1};
    if (length(varargin) < 2)
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{2};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    analysis = viewGet(view,'analysis',analysisNum);
    % handle if passed in a number
    if isnumeric(overlayName)
      if any((overlayName > 0) & (overlayName <= viewGet(view,'numOverlays')))
	val = overlayName(logical(overlayName));
      end
    elseif ~isempty(analysis) & ~isempty(analysis.overlays)
      overlayNames = {analysis.overlays(:).name};
      val = find(ismember(overlayNames,overlayName));
    end
  case {'overlay'}
    % overlay = viewGet(view,'overlay',[overlayNum],[analysisNum])
    % overlay = viewGet(view,'overlay',overlayNum,[])
    % overlay = viewGet(view,'overlay',[],analysisNum)
    % overlay = viewGet(view,'overlay',[],[])
    % overlay = viewGet(view,'overlay',overlayNum)
    % overlay = viewGet(view,'overlay')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    switch (length(varargin))
      case 1
        overlayNum = varargin{1};
        analysisNum = viewGet(view,'currentAnalysis');
      case 2
        overlayNum = varargin{1};
        analysisNum = varargin{2};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    if isempty(overlayNum)
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    % if overlayNum is a string, then convert to a number
    if isstr(overlayNum),overlayNum = viewGet(view,'overlayNum',overlayNum);end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis) & ~isempty(analysis.overlays)
      n = viewGet(view,'numberofOverlays',analysisNum);
      if overlayNum & (overlayNum > 0) & (overlayNum <= n)
        val = analysis.overlays(overlayNum);
      end
    end
    if isempty(val),disp(sprintf('(viewGet:overlay) Could not find requested overlay.'));end
  case {'overlaynames'}
    % overlayNames = viewGet(view,'overlayNames',[analysisNum])
    % overlayNames = viewGet(view,'overlayNames',[])
    % overlayNames = viewGet(view,'overlayNames')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
    else
      analysisNum = varargin{1};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(view.analyses) & ~isempty(analysis.overlays)
      val = {view.analyses{analysisNum}.overlays(:).name};
    end
  case {'overlayname'}
    % overlayName = viewGet(view,'overlayName',[overlayNum],[analysisNum])
    % overlayName = viewGet(view,'overlayName',overlayNum,[])
    % overlayName = viewGet(view,'overlayName',[],analysisNum)
    % overlayName = viewGet(view,'overlayName',[],[])
    % overlayName = viewGet(view,'overlayName',overlayNum)
    % overlayName = viewGet(view,'overlayName')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    switch (length(varargin))
      case 1
        overlayNum = varargin{1};
        analysisNum = viewGet(view,'currentAnalysis');
      case 2
        overlayNum = varargin{1};
        analysisNum = varargin{2};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    if isempty(overlayNum)
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if isstr(overlayNum),overlayNum = viewGet(view,'overlayNum',overlayNum);end
    if ~isempty(analysis) & ~isempty(analysis.overlays)
      n = viewGet(view,'numberofOverlays',analysisNum);
      if overlayNum & (overlayNum > 0) & (overlayNum <= n)
        val = analysis.overlays(overlayNum).name;
      end
    end
  case {'overlaytype'}
    % overlayType = viewGet(view,'overlayType',[overlayNum],[analysisNum])
    % overlayType = viewGet(view,'overlayType',overlayNum,[])
    % overlayType = viewGet(view,'overlayType',[],analysisNum)
    % overlayType = viewGet(view,'overlayType',[],[])
    % overlayType = viewGet(view,'overlayType',overlayNum)
    % overlayType = viewGet(view,'overlayType')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    switch (length(varargin))
      case 1
        overlayNum = varargin{1};
        analysisNum = viewGet(view,'currentAnalysis');
      case 2
        overlayNum = varargin{1};
        analysisNum = varargin{2};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    if isempty(overlayNum)
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis) & ~isempty(analysis.overlays)
      n = viewGet(view,'numberofOverlays',analysisNum);
      if overlayNum & (overlayNum > 0) & (overlayNum <= n)
        val = analysis.overlays(overlayNum).type;
      end
    end
  case{'overlayxform'}
    % xform = viewGet(view,'overlayXform',scanNum)
    if ieNotDefined('varargin')
      mrErrorDlg('viewGet overlayXform: must specify scan.');
    end
    scan = varargin{1};
    val = viewGet(view,'scanxform',scan);
  case {'overlaydataval'}
    % overlaydata = viewGet(view,'overlayDataVal',x,y,s)
    % gets the value of the current overlay at the x,y,s point
    val = [];
    if ~isempty(view.curAnalysis)
      curOverlay = view.analyses{view.curAnalysis}.curOverlay;
      if ~isempty(curOverlay)
	overlay = view.analyses{view.curAnalysis}.overlays(curOverlay).data;
	if ~isempty(overlay{view.curScan})
	  val = overlay{view.curScan}(varargin{1},varargin{2},varargin{3});
	end
      end
    end
  case {'overlaydata'}
    % overlaydata = viewGet(view,'overlaydata',scanNum,[overlayNum],[analysisNum])
    % overlaydata = viewGet(view,'overlaydata',scanNum,overlayNum,[])
    % overlaydata = viewGet(view,'overlaydata',scanNum,[],analysisNum)
    % overlaydata = viewGet(view,'overlaydata',scanNum,[],[])
    % overlaydata = viewGet(view,'overlaydata',scanNum,overlayNum)
    % overlaydata = viewGet(view,'overlaydata',scanNum)
    % scanNum can be a number or []. If [] then it returns the entire cell
    % array of data.
    if ieNotDefined('varargin')
      mrErrorDlg('viewGet overlayData: must specify scan.');
    end
    scan = varargin{1};
    switch (length(varargin))
      case 1
        analysisNum = viewGet(view,'currentAnalysis');
        overlayNum = viewGet(view,'currentOverlay',analysisNum);
      case 2
        overlayNum = varargin{2};
        analysisNum = viewGet(view,'currentAnalysis');
      case 3
        overlayNum = varargin{2};
        analysisNum = varargin{3};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    if isempty(overlayNum)
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis) & ~isempty(analysis.overlays)
      n = viewGet(view,'numberofOverlays',analysisNum);
      if all(overlayNum & (overlayNum > 0) & (overlayNum <= n))
        if isempty(scan) %if we want the whole data structure
          scanList = 1:length(analysis.overlays(overlayNum(end)).data); %we assume that the number of scans is identical for all overlays
        else
          scanList = scan;
        end
        cScan = 0;
        for iScan = scanList
          cScan = cScan+1;
          scanDims = viewGet(view,'scanDims',iScan);
          cOverlay=0;
          for iOverlay = overlayNum
            cOverlay = cOverlay+1;
            data = analysis.overlays(iOverlay).data;
            if length(data) >= iScan && ~isempty(data{iScan})
              val{cScan}(:,:,:,cOverlay) = data{iScan};
            else
	      %  No overlay found, filling with nans. check for empty scan dims
	      % if no scandims return empty data.
	      if isempty(scanDims)
		val{cScan} = [];
	      else
		val{cScan}(:,:,:,cOverlay) = NaN(scanDims);
	      end
            end
          end
        end
        if ~isempty(scan)
          val = val{1};
        end
      end
    end
    
  case  {'maxoverlaydata'}
    % maxoverlaydata = viewGet(view,'maxoverlaydata',[overlayNum],[analysisNum],[scanlist])
    switch (length(varargin))
      case 0
        analysisNum = viewGet(view,'currentAnalysis');
        overlayNum = viewGet(view,'currentOverlay');
      case 1
        overlayNum = varargin{1};
        analysisNum = viewGet(view,'currentAnalysis');
      case {2,3}
        overlayNum = varargin{1};
        analysisNum = varargin{2};
    end
    overlay = viewGet(view,'overlay',overlayNum,analysisNum);
% %     if isfield(overlay,'maxOverlayData')     % This was meant to avoid having to recompute minoverlaydata that has been compute before
% %       val = overlay.maxoverlaydata;          % but turns out to be too much of a headache, plus what if an overlay is added to a scan and the value changes ?
% %     else
    if ~isempty(overlay)
      switch (length(varargin))
        case {0 1 2}
          scanlist = 1:length(overlay.data);
        case 3
          scanlist = varargin{3};
          scanlist = scanlist(scanlist<=length(overlay.data));
      end
      val = -inf;
      for iOverlay = scanlist
         if ~isempty(overlay.data{iOverlay})
            val = max(val,max(overlay.data{iOverlay}(:)));
         end
      end
      if val==-inf
         val = [];
      end
    end
% %       viewSet(view,'maxoverlaydata',val,overlayNum);
% %     end
    
  case  {'minoverlaydata'}
    % minoverlaydata = viewGet(view,'minoverlaydata',[overlayNum],[analysisNum],[scanlist])
    switch (length(varargin))
      case 0
        analysisNum = viewGet(view,'currentAnalysis');
        overlayNum = viewGet(view,'currentOverlay');
      case 1
        overlayNum = varargin{1};
        analysisNum = viewGet(view,'currentAnalysis');
      case {2,3}
        overlayNum = varargin{1};
        analysisNum = varargin{2};
    end
    overlay = viewGet(view,'overlay',overlayNum,analysisNum);
% %     if isfield(overlay,'minoverlaydata')    % This was meant to avoid having to recompute minoverlaydata that has been compute before
% %       val = overlay.minoverlaydata;         % but turns out to be too much of a headache, plus what if an overlay is added to a scan and the value changes ?
% %     else
    if ~isempty(overlay)
      switch (length(varargin))
        case {0 1 2}
          scanlist = 1:length(overlay.data);
        case 3
          scanlist = varargin{3};
          scanlist = scanlist(scanlist<=length(overlay.data));
      end
      val = inf;
      for iOverlay = scanlist
         if ~isempty(overlay.data{iOverlay})
            val = min(val,min(overlay.data{iOverlay}(:)));
         end
      end
      if val==inf
         val = [];
      end
    end
% %       viewSet(view,'minoverlaydata',val,overlayNum);
% %     end
    
  case {'overlaydims'}
    % overlaydims = viewGet(view,'overlaydims',scanNum,[overlayNum],[analysisNum])
    % overlaydims = viewGet(view,'overlaydims',scanNum,overlayNum,[])
    % overlaydims = viewGet(view,'overlaydims',scanNum,[],analysisNum)
    % overlaydims = viewGet(view,'overlaydims',scanNum,[],[])
    % overlaydims = viewGet(view,'overlaydims',scanNum,overlayNum)
    % overlaydims = viewGet(view,'overlaydims',scanNum)
    if ieNotDefined('varargin')
      mrErrorDlg('viewGet overlayDims: must specify scan.');
    end
    scan = varargin{1};
    switch (length(varargin))
      case 1
        analysisNum = viewGet(view,'currentAnalysis');
        overlayNum = viewGet(view,'currentOverlay',analysisNum);
      case 2
        overlayNum = varargin{2};
        analysisNum = viewGet(view,'currentAnalysis');
      case 3
        overlayNum = varargin{2};
        analysisNum = varargin{3};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    if isempty(overlayNum)
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis) & ~isempty(analysis.overlays)
      n = viewGet(view,'numberofOverlays',analysisNum);
      if overlayNum & (overlayNum > 0) & (overlayNum <= n)
        val = size(analysis.overlays(overlayNum).data{scan});
      end
    end
  case {'overlayclip'}
    % overlayclip = viewGet(view,'overlayclip',[overlayNum],[analysisNum])
    % overlayclip = viewGet(view,'overlayclip',overlayNum,[])
    % overlayclip = viewGet(view,'overlayclip',[],analysisNum)
    % overlayclip = viewGet(view,'overlayclip',[],[])
    % overlayclip = viewGet(view,'overlayclip',overlayNum)
    % overlayclip = viewGet(view,'overlayclip')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    switch (length(varargin))
      case 1
        overlayNum = varargin{1};
        analysisNum = viewGet(view,'currentAnalysis');
      case 2
        overlayNum = varargin{1};
        analysisNum = varargin{2};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    if isempty(overlayNum)
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis) & ~isempty(analysis.overlays)
      n = viewGet(view,'numberofOverlays',analysisNum);
      if overlayNum & (overlayNum > 0) & (overlayNum <= n)
        val = analysis.overlays(overlayNum).clip;
      end
    end
  case {'alphaoverlay'}
    % alphaoverlay = viewGet(view,'alphaOverlay',[overlaynum])
    analysisNum = viewGet(view,'currentAnalysis');
    if isempty(varargin) || isempty(varargin{1})
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    else
      overlayNum = varargin{1};
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis) & ~isempty(analysis.overlays)
      n = viewGet(view,'numberofOverlays',analysisNum);
      if overlayNum & (overlayNum > 0) & (overlayNum <= n)
        for iOverlay=1:length(overlayNum)
          val{iOverlay} = char(analysis.overlays(overlayNum(iOverlay)).alphaOverlay);
        end
      end
      if length(val)==1
        val = val{1};
      end
    end
  case {'alphaoverlayexponent'}
    % overlayclip = viewGet(view,'alphaOverlayExponent',[overlaynum])
    val = 1;
    analysisNum = viewGet(view,'currentAnalysis');
    if isempty(varargin) || isempty(varargin{1})
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    else
      overlayNum = varargin{1};
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis) & ~isempty(analysis.overlays)
      n = viewGet(view,'numberofOverlays',analysisNum);
      if overlayNum & all(overlayNum > 0) & all(overlayNum <= n)
        val = [analysis.overlays(overlayNum).alphaOverlayExponent];
      end
    end
  case {'overlaycmap'}
    % overlaycmap = viewGet(view,'overlaycmap',[overlayNum],[analysisNum])
    % overlaycmap = viewGet(view,'overlaycmap',overlayNum,[])
    % overlaycmap = viewGet(view,'overlaycmap',[],analysisNum)
    % overlaycmap = viewGet(view,'overlaycmap',[],[])
    % overlaycmap = viewGet(view,'overlaycmap',overlayNum)
    % overlaycmap = viewGet(view,'overlaycmap')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    switch (length(varargin))
      case 1
        overlayNum = varargin{1};
        analysisNum = viewGet(view,'currentAnalysis');
      case 2
        overlayNum = varargin{1};
        analysisNum = varargin{2};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    if isempty(overlayNum)
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis) & ~isempty(analysis.overlays)
      n = viewGet(view,'numberofOverlays',analysisNum);
      if overlayNum & (overlayNum > 0) & (overlayNum <= n)
        val = analysis.overlays(overlayNum).colormap;
      end
    end
  case {'overlayctype'}
    % overlayctype = viewGet(view,'overlayctype',[overlayNum],[analysisNum])
    % overlayctype = viewGet(view,'overlayctype',overlayNum,[])
    % overlayctype = viewGet(view,'overlayctype',[],analysisNum)
    % overlayctype = viewGet(view,'overlayctype',[],[])
    % overlayctype = viewGet(view,'overlayctype',overlayNum)
    % overlayctype = viewGet(view,'overlayctype')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    switch (length(varargin))
      case 1
        overlayNum = varargin{1};
        analysisNum = viewGet(view,'currentAnalysis');
      case 2
        overlayNum = varargin{1};
        analysisNum = varargin{2};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    if isempty(overlayNum)
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis) & ~isempty(analysis.overlays)
      n = viewGet(view,'numberofOverlays',analysisNum);
      if overlayNum & (overlayNum > 0) & (overlayNum <= n)
        if isfield(analysis.overlays(overlayNum),'colormapType')
          val = analysis.overlays(overlayNum).colormapType;
        else
          val = 'normal';
        end
      end
    end
  case {'overlayrange'}
    % overlayrange = viewGet(view,'overlayrange',[overlayNum],[analysisNum])
    % overlayrange = viewGet(view,'overlayrange',overlayNum,[])
    % overlayrange = viewGet(view,'overlayrange',[],analysisNum)
    % overlayrange = viewGet(view,'overlayrange',[],[])
    % overlayrange = viewGet(view,'overlayrange',overlayNum)
    % overlayrange = viewGet(view,'overlayrange')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    switch (length(varargin))
      case 1
        overlayNum = varargin{1};
        analysisNum = viewGet(view,'currentAnalysis');
      case 2
        overlayNum = varargin{1};
        analysisNum = varargin{2};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    if isempty(overlayNum)
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis) & ~isempty(analysis.overlays)
      n = viewGet(view,'numberofOverlays',analysisNum);
      if overlayNum & (overlayNum > 0) & (overlayNum <= n)
        val = analysis.overlays(overlayNum).range;
      end
    end
    
  case {'overlayalpha'}
    % overlayalpha = viewGet(view,'overlayalpha',[overlayNum],[analysisNum])
    % overlayalpha = viewGet(view,'overlayalpha',overlayNum,[])
    % overlayalpha = viewGet(view,'overlayalpha',[],analysisNum)
    % overlayalpha = viewGet(view,'overlayalpha',[],[])
    % overlayalpha = viewGet(view,'overlayalpha',overlayNum)
    % overlayalpha = viewGet(view,'overlayalpha')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    switch (length(varargin))
      case 1
        overlayNum = varargin{1};
        analysisNum = viewGet(view,'currentAnalysis');
      case 2
        overlayNum = varargin{1};
        analysisNum = varargin{2};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    if isempty(overlayNum)
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis) & ~isempty(analysis.overlays)
      n = viewGet(view,'numberofOverlays',analysisNum);
      if overlayNum & (overlayNum > 0) & (overlayNum <= n)
        val = analysis.overlays(overlayNum).alpha;
      end
    end
  case {'overlaydate'}
    % overlaydate = viewGet(view,'overlaydate',[overlayNum],[analysisNum])
    % overlaydate = viewGet(view,'overlaydate',overlayNum,[])
    % overlaydate = viewGet(view,'overlaydate',[],analysisNum)
    % overlaydate = viewGet(view,'overlaydate',[],[])
    % overlaydate = viewGet(view,'overlaydate',overlayNum)
    % overlaydate = viewGet(view,'overlaydate')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    switch (length(varargin))
      case 1
        overlayNum = varargin{1};
        analysisNum = viewGet(view,'currentAnalysis');
      case 2
        overlayNum = varargin{1};
        analysisNum = varargin{2};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    if isempty(overlayNum)
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis) & ~isempty(analysis.overlays)
      n = viewGet(view,'numberofOverlays',analysisNum);
      if overlayNum & (overlayNum > 0) & (overlayNum <= n)
        val = analysis.overlays(overlayNum).date;
      end
    end
  case {'interrogator'}
    % interrogator = viewGet(view,'interrogator',[overlayNum],[analysisNum])
    % interrogator = viewGet(view,'interrogator',overlayNum,[])
    % interrogator = viewGet(view,'interrogator',[],analysisNum)
    % interrogator = viewGet(view,'interrogator',[],[])
    % interrogator = viewGet(view,'interrogator',overlayNum)
    % interrogator = viewGet(view,'interrogator')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    switch (length(varargin))
      case 1
        overlayNum = varargin{1};
        analysisNum = viewGet(view,'currentAnalysis');
      case 2
        overlayNum = varargin{1};
        analysisNum = varargin{2};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    if isempty(overlayNum)
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis) & ~isempty(analysis.overlays)
      n = viewGet(view,'numberofOverlays',analysisNum);
      if overlayNum & (overlayNum > 0) & (overlayNum <= n)
        val = analysis.overlays(overlayNum).interrogator;
      end
    end
    if isempty(val)
      val = 'timecoursePlot';
    end
  case {'overlaygroupname'}
    % groupName = viewGet(view,'overlayGroupName',[overlayNum],[analysisNum])
    % groupName = viewGet(view,'overlayGroupName',overlayNum,[])
    % groupName = viewGet(view,'overlayGroupName',[],analysisNum)
    % groupName = viewGet(view,'overlayGroupName',[],[])
    % groupName = viewGet(view,'overlayGroupName',overlayNum)
    % groupName = viewGet(view,'overlayGroupName')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    switch (length(varargin))
      case 1
        overlayNum = varargin{1};
        analysisNum = viewGet(view,'currentAnalysis');
      case 2
        overlayNum = varargin{1};
        analysisNum = varargin{2};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    if isempty(overlayNum)
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis) & ~isempty(analysis.overlays)
      n = viewGet(view,'numberofOverlays',analysisNum);
      if overlayNum & (overlayNum > 0) & (overlayNum <= n)
        val = analysis.overlays(overlayNum).groupName;
      end
    end
  case {'overlayparams'}
    % params = viewGet(view,'overlayParams',[overlayNum],[analysisNum])
    % params = viewGet(view,'overlayParams',overlayNum,[])
    % params = viewGet(view,'overlayParams',[],analysisNum)
    % params = viewGet(view,'overlayParams',[],[])
    % params = viewGet(view,'overlayParams',overlayNum)
    % params = viewGet(view,'overlayParams')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    switch (length(varargin))
      case 1
        overlayNum = varargin{1};
        analysisNum = viewGet(view,'currentAnalysis');
      case 2
        overlayNum = varargin{1};
        analysisNum = varargin{2};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    if isempty(overlayNum)
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis) & ~isempty(analysis.overlays)
      n = viewGet(view,'numberofOverlays',analysisNum);
      if overlayNum & (overlayNum > 0) & (overlayNum <= n)
        val = analysis.overlays(overlayNum).params;
      end
    end
  case {'overlayreconcilefunction'}
    % reconcileFunction = viewGet(view,'overlayReconcileFunction',[overlayNum],[analysisNum])
    % reconcileFunction = viewGet(view,'overlayReconcileFunction',overlayNum,[])
    % reconcileFunction = viewGet(view,'overlayReconcileFunction',[],analysisNum)
    % reconcileFunction = viewGet(view,'overlayReconcileFunction',[],[])
    % reconcileFunction = viewGet(view,'overlayReconcileFunction',overlayNum)
    % reconcileFunction = viewGet(view,'overlayReconcileFunction')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    switch (length(varargin))
      case 1
        overlayNum = varargin{1};
        analysisNum = viewGet(view,'currentAnalysis');
      case 2
        overlayNum = varargin{1};
        analysisNum = varargin{2};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    if isempty(overlayNum)
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis) & ~isempty(analysis.overlays)
      n = viewGet(view,'numberofOverlays',analysisNum);
      if overlayNum & (overlayNum > 0) & (overlayNum <= n)
        val = analysis.overlays(overlayNum).reconcileFunction;
      end
    end
  case {'overlaymergefunction'}
    % mergeFunction = viewGet(view,'overlayMergeFunction',[overlayNum],[analysisNum])
    % mergeFunction = viewGet(view,'overlayMergeFunction',overlayNum,[])
    % mergeFunction = viewGet(view,'overlayMergeFunction',[],analysisNum)
    % mergeFunction = viewGet(view,'overlayMergeFunction',[],[])
    % mergeFunction = viewGet(view,'overlayMergeFunction',overlayNum)
    % mergeFunction = viewGet(view,'overlayMergeFunction')
    if ieNotDefined('varargin')
      analysisNum = viewGet(view,'currentAnalysis');
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    switch (length(varargin))
      case 1
        overlayNum = varargin{1};
        analysisNum = viewGet(view,'currentAnalysis');
      case 2
        overlayNum = varargin{1};
        analysisNum = varargin{2};
    end
    if isempty(analysisNum)
      analysisNum = viewGet(view,'currentAnalysis');
    end
    if isempty(overlayNum)
      overlayNum = viewGet(view,'currentOverlay',analysisNum);
    end
    analysis = viewGet(view,'analysis',analysisNum);
    if ~isempty(analysis) & ~isempty(analysis.overlays)
      n = viewGet(view,'numberofOverlays',analysisNum);
      if overlayNum & (overlayNum > 0) & (overlayNum <= n)
        val = analysis.overlays(overlayNum).mergeFunction;
      end
    end
    
    % corAnal
  case {'coranal'}
    % corAnal = viewGet(view,'corAnal')
    % Find all corAnals
    analysisTypes = viewGet(view,'analysisTypes');
    corAnalNums = find(strcmp('corAnal',analysisTypes));
    if ~isempty(corAnalNums)
      % Check if current analysis is corAnal. Otherwise, use the
      % first corAnal on the analyses list.
      curAnalysisNum = viewGet(view,'currentAnalysis');
      if find(curAnalysisNum == corAnalNums)
        n = curAnalysisNum;
      else
        n = corAnalNums(1);
      end
      val = viewGet(view,'analysis',n);
    end
  case {'coranalparams'}
    % corAnalParams = viewGet(view,'corAnalParams')
    corAnal = viewGet(view,'corAnal');
    if ~isempty(corAnal)
      val = corAnal.params;
    end
  case {'co'}
    % co = viewGet(view,'co')
    n = viewGet(view,'overlayNum','co');
    if n
      val = viewGet(view,'overlay',n);
    end
  case {'amp'}
    % amp = viewGet(view,'amp')
    n = viewGet(view,'overlayNum','amp');
    if n
      val = viewGet(view,'overlay',n);
    end
  case {'ph'}
    % ph = viewGet(view,'ph')
    n = viewGet(view,'overlayNum','ph');
    if n
      val = viewGet(view,'overlay',n);
    end
  case {'ncycles'}
    % ncycles = viewGet(view,'ncycles',scanNum)
    if ieNotDefined('varargin')
      mrErrorDlg('viewGet ncycles: must specify scan.');
    end
    scanNum = varargin{1};
    params = viewGet(view,'corAnalParams');
    if ~isempty(params)
      if isfield(params,'ncycles')
        val = params.ncycles(scanNum);
      else
        val = params.scanParams{scanNum}.ncycles;
      end
    end
  case {'detrend'}
    % detrend = viewGet(view,'detrend',scanNum)
    if ieNotDefined('varargin')
      mrErrorDlg('viewGet detrend: must specify scan.');
    end
    scanNum = varargin{1};
    params = viewGet(view,'corAnalParams');
    if ~isempty(params)
      if isfield(params,'detrend')
        val = params.detrend{scanNum};
      else
        val = params.scanParams{scanNum}.detrend;
      end
    end
  case {'spatialnorm'}
    % spatialnorm = viewGet(view,'spatialnorm',scanNum)
    if ieNotDefined('varargin')
      mrErrorDlg('viewGet spatialnorm: must specify scan.');
    end
    scanNum = varargin{1};
    params = viewGet(view,'corAnalParams');
    if ~isempty(params)
      if isfield(params,'spatialnorm')
        val = params.spatialnorm{scanNum};
      else
        val = params.scanParams{scanNum}.spatialnorm;
      end
    end
  case {'trigonometricfunction'}
    % spatialnorm = viewGet(view,'spatialnorm',scanNum)
    if ieNotDefined('varargin')
      mrErrorDlg('viewGet trigonometricFunction: must specify scan.');
    end
    scanNum = varargin{1};
    params = viewGet(view,'corAnalParams');
    if ~isempty(params)
      if isfield(params,'trigonometricFunction')
        val = params.trigonometricFunction{scanNum};
      elseif isfield(params,'scanParams') && isfield(params.scanParams{scanNum},'trigonometricFunction')
        val = params.scanParams{scanNum}.trigonometricFunction;
      else
        val = 'Sine'; %value by default because it used to not be a parameter
      end
    end

    
    % GUI
  case {'fignum','figurenumber'}
    % figureHandle = viewGet(view,'figureNumber');
    val = view.figure;
  case {'curscan','currentscan'}
    % scan = viewGet(view,'currentScan');
    val = view.curScan;
    if isempty(val),val = 1;end
  case {'curslice','currentslice'}
    % slice = viewGet(view,'currentSlice');
    if isfield(view.curslice,'sliceNum')
      val = view.curslice.sliceNum;
    elseif viewGet(view,'baseMultiAxis')==0
      % default to 1
      val = 1;
    end
  case {'curcoords','currentcoordinates'}
    % coords = viewGet(view,'currentCoordinates');
    fig = viewGet(view,'fignum');
    if ~isempty(fig)
      handles = guidata(fig);
      val = handles.coords;
      % never let coords be 0
      val(val==0) = 1;
    else
      val = [];
    end
  case {'alpha'}
    % alpha = viewGet(view,'alpha',[overlaynum]);
    if isempty(varargin) || isempty(varargin{1})
      fig = viewGet(view,'fignum');
      overlayNum = viewGet(view,'currentOverlay');
    else
      fig =[];
      overlayNum = varargin{1};
    end
    if ~isempty(fig)
      handles = guidata(fig);
      val = get(handles.alphaSlider,'Value');
    else
      % get alpha from analysis structure
      analysisNum = viewGet(view,'currentAnalysis');
      if ~isempty(analysisNum) & ~isempty(overlayNum) &  ~isempty(view.analyses{analysisNum}.overlays)
        val = view.analyses{analysisNum}.overlays(overlayNum).alpha;
      else
        val = 1;
      end
    end
  case {'overlaymin'}
    % overlayMin = viewGet(view,'overlayMin');
    % overlayMin = viewGet(view,'overlayMin',<overlayName>);
    fig = viewGet(view,'fignum'); %JB: this is a quick fix for the case in which several overlays are loaded
    if ~isempty(fig)              %    but it doesn't solve the problem if used form a script (no figure)
      handles = guidata(fig);          
      val = get(handles.overlayMinSlider,'Value');
    else
      if (length(varargin) > 0)
        overlayNum = viewGet(view,'overlayNum',varargin{1});
      else
        overlayNum = viewGet(view,'currentOverlay');
      end
      if isempty(overlayNum)
        disp(sprintf('(viewGet:overlayMin) Unknown overlay'));
      end
      % get overlayMin from analysis structure
      analysisNum = viewGet(view,'currentAnalysis');
      if ~isempty(analysisNum) & ~isempty(overlayNum) &  ~isempty(view.analyses{analysisNum}.overlays)
        val = view.analyses{analysisNum}.overlays(overlayNum).clip(1);
      else
        val = 0;
      end
    end
  case {'overlaymax'}
    % overlayMax = viewGet(view,'overlayMax');
    % overlayMax = viewGet(view,'overlayMax',<overlayName>);
    fig = viewGet(view,'fignum'); %JB: this is a quick fix for the case in which several overlays are loaded
    if ~isempty(fig)              %    but it doesn't solve the problem if used from a script (no figure)
      handles = guidata(fig);
      val = get(handles.overlayMaxSlider,'Value');
    else
      if (length(varargin) > 0)
        overlayNum = viewGet(view,'overlayNum',varargin{1});
      else
        overlayNum = viewGet(view,'currentOverlay');
      end
      if isempty(overlayNum)
        disp(sprintf('(viewGet:overlayMin) Unknown overlay'));
      end
      % get overlayMax from analysis structure
      analysisNum = viewGet(view,'currentAnalysis');
      if ~isempty(analysisNum) & ~isempty(overlayNum) &  ~isempty(view.analyses{analysisNum}.overlays)
        val = view.analyses{analysisNum}.overlays(overlayNum).clip(2);
      else
        val = 1;
      end
    end
  case {'rotate'}
    % rotate = viewGet(view,'rotate');
    baseType = viewGet(view,'baseType');
    % rotation is dependent on the base type
    % since surfaces aren't rotated in the
    % same way, they are rotated through
    % rotateSurface
    if baseType <= 1
        curBase = viewGet(view,'curBase');
        if ~isempty(curBase)
          val = view.baseVolumes(curBase).rotate;
        else
          val = 0;
        end
    else
      val = 0;
    end
  case {'rotatesurface'}
    % rotate = viewGet(view,'rotateSurface');
    baseType = viewGet(view,'baseType');
    if baseType == 2
      fig = viewGet(view,'fignum');
      if ~isempty(fig)
        handles = guidata(fig);
        % 360- is so that the anatomy rotates correct direction
        val = 360-get(handles.rotateSlider,'Value');
      else
        val = 0;
      end
    elseif baseType == 0
      val = 0;
      % normally does not apply to inplanes, but if we are displaying
      % all planes of anatomy then it applies to the 3D slice view
      if viewGet(view,'baseMultiAxis')>0
	fig = viewGet(view,'fignum');
	if ~isempty(fig)
	  handles = guidata(fig);
	  val = 360-get(handles.rotateSlider,'Value');
	end
      end	  
    end
    
  case {'corticaldepthbins'}
    % corticaldepthBins = viewGet(view,'corticaldepthBins');
    val = mrGetPref('corticalDepthBins');
    
  case {'corticaldepth'}
    % corticaldepth = viewGet(view,'corticaldepth');
    val = sort([viewGet(view,'corticalMinDepth') viewGet(view,'corticalMaxDepth')]);
    
  case {'corticalmindepth'}
    % corticalmindepth = viewGet(view,'corticalmindepth');
    if isfield(view.curslice,'corticalDepth')
      val = min(view.curslice.corticalDepth);
    else
      fig = viewGet(view,'fignum');
      if ~isempty(fig)
        handles = guidata(fig);
        val = get(handles.corticalDepthSlider,'Value');
    else
      corticalDepthBins=viewGet(view,'corticaldepthbins');
      val=round((corticalDepthBins-1)/2)/(corticalDepthBins-1);
      end
    end
    
  case {'corticalmaxdepth'}
    % corticalmaxdepth = viewGet(view,'corticalmaxdepth');
    if isfield(view.curslice,'corticalDepth')
      val = max(view.curslice.corticalDepth);
    else
      fig = viewGet(view,'fignum');
      if ~isempty(fig)
        handles = guidata(fig);
        if isfield(handles,'corticalMaxDepthSlider')
          val = get(handles.corticalMaxDepthSlider,'Value');
        else
          val = [];
        end
      else
        val = [];
      end
    end
    
  case {'sliceorientation'}
    % sliceorientation = viewGet(view,'sliceorientation');
    val = view.sliceOrientation;
  case {'cursliceoverlaycoords'}
    % overlayCoords = viewGet(view,'cursliceoverlaycoords');
    val = view.curslice.overlayCoords;
  case {'curslicebasecoords'}
    % baseCoords = viewGet(view,'curslicebasecoords');
    val = view.curslice.baseCoords;
  case {'defaultinterrogators'}
    % defaultInterrogators = viewGet(view,'defaultInterrogators');
    % see the viewSet for more info
    if isfield(MLR,'defaultInterrogators')
      val = MLR.defaultInterrogators;
    else
      val = [];
    end
  case {'colormaps'}
    % colormaps = viewGet(view,'colormaps');
    % see the viewSet for more info
    if isfield(MLR,'colormaps')
      val = MLR.colormaps;
    else
      val = [];
    end
    
  case{'clipacrossoverlays'}
    % clipacrossoverlays = viewGet(view,'clipacrossoverlays',[analysisNum]);
    if (length(varargin) > 0)
      analysisNum = varargin{1};
    else
      analysisNum = viewGet(view,'currentAnalysis');
    end
    val = true;
    if ~isempty(analysisNum) && (analysisNum > 0) && (analysisNum <= length(view.analyses))
      val = view.analyses{analysisNum}.clipAcrossOverlays;
    end

  case{'colorblending'}
    val='Alpha blend';
    if isfield(view,'figure')
      handles = guidata(view.figure);
      if isfield(handles,'colorBlendingPopup') 
        colorBlendingOptions = get(handles.colorBlendingPopup,'string');
        val = colorBlendingOptions{get(handles.colorBlendingPopup,'value')};
      end
    end

    
 otherwise
    if isempty(view)
      disp(sprintf('(viewGet) No viewGet for %s',param));
      return
    else
      switch(lower(view.viewType))
        case 'volume'
          val = volumeGet(view,param,varargin{:});
        otherwise
          error('Unknown type of View.');
      end
    end
end

return;


% N.B.  For all the routines below, varargin{1} is the varargin passed in by the
% main routine. Hence, the first entry varargin{1}, is the entire varargin
% in the calling routine.
%------------------------------

%------------------------------
function val  = volumeGet(view,param,varargin)

global MLR

val = [];
switch lower(param)
  
  case {'tseriessize'}
    % size of entire tseries: [ysize xsize zsize nframes]
    %
    % tseriessize = viewGet(view,'tseriessize',[scanNum],[groupNum])
    % tseriessize = viewGet(view,'tseriessize',scanNum,[])
    % tseriessize = viewGet(view,'tseriessize',scanNum)
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      datasize = viewGet(view,'datasize',s,g);
      nframes = viewGet(view,'nframes',s,g);
      val = [datasize nframes];
    end
    
  case {'datasize','dims'}
    % dims of single temporal frame of functional volume (same as size
    % of parameter map) for a given scan.
    %
    % datasize = viewGet(view,'datasize',[scanNum],[groupNum])
    % datasize = viewGet(view,'datasize',scanNum,[])
    % datasize = viewGet(view,'datasize',scanNum)
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      val = MLR.groups(g).scanParams(s).dataSize;
    end
    
  case {'slicedims'}
    % dims of single slice and single temporal frame of functional
    % volume (same as dims of single slice of parameter map) for a
    % given scan.
    %
    % slicedims = viewGet(view,'slicedims',[scanNum],[groupNum])
    % slicedims = viewGet(view,'slicedims',scanNum,[])
    % slicedims = viewGet(view,'slicedims',scanNum)
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      datasize = viewGet(view,'datasize',s,g);
      val = datasize(1:2);
    end
    
  case {'nslices'}
    % nslices = viewGet(view,'nslices',[scanNum],[groupNum])
    % nslices = viewGet(view,'nslices',scanNum,[])
    % nslices = viewGet(view,'nslices',scanNum)
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      datasize = viewGet(view,'datasize',s,g);
      val = datasize(3);
    end
    
  case{'sliceorder'}
    % n = viewGet(view,'sliceOrder',[scanNum],[groupNum])
    % n = viewGet(view,'sliceOrder',scanNum,[])
    % n = viewGet(view,'sliceOrder',scanNum)
    [s g] = getScanAndGroup(view,varargin,param);
    nscans = viewGet(view,'nscans',g);
    if (nscans >= s) & (s > 0)
      nslices = viewGet(view,'nslices',s,g);
      if strcmp(mrGetPref('site'),'NYU')
        if isodd(nslices)
          val = [1:2:nslices,2:2:nslices-1];
        else
          val = [2:2:nslices,1:2:nslices-1];
        end
      elseif strcmp(lower(mrGetPref('site')),'nottingham')
        %This is for the interleave setting on a Philips scanner (7T)
        jump=round(sqrt(nslices));  %the jump is the closest integer to the square root of the number of slices
        val=[];
        for i=1:jump
          val = [val i:jump:nslices];
        end
      else
	% check for fidinfo
	fidInfo = viewGet(view,'fidInfo',s,g);
	if isempty(fidInfo) 
	  % DEFAULT, warn user and return slices in slice order
	  mrWarnDlg('(viewGet) Slice ordering is unknown for this site. Using default order: [1:nslices]. If this is incorrect, then edit viewGet sliceOrder to add the convention for your site.');
	  val = [1:nslices];
	else
	  % extract from fidInfo
	  if length(fidInfo) > 1
	    disp(sprintf('(viewGet:sliceOrder) There seems to be more than one associated FIDs with this scan, returning the sliceOrder for the first fid: %s',fidInfo{1}.fidname));
	  end
	  % get slice order
	  val = fidInfo{1}.sliceOrder;
	end
      end
    end
    
  case{'slicetimes'}
    % n = viewGet(view,'sliceTimes',[scanNum],[groupNum])
    % n = viewGet(view,'sliceTimes',scanNum,[])
    % n = viewGet(view,'sliceTimes',scanNum)
    [s g] = getScanAndGroup(view,varargin,param);
    sliceOrder = viewGet(view,'sliceOrder',s,g);
    nslices = viewGet(view,'nslices',s,g);
    val = zeros(1,nslices);
    val(sliceOrder) = [0:nslices-1]/nslices;
  case{'panelhandle'}
   % h = viewGet(view,'panelHandle','Panel name')
   val = [];
   if length(varargin)<1
     disp(sprintf('(viewGet) Must specify panel name to get panel handle'));
     return
   end
   if isfield(MLR,'panels')
     for iPanel = 1:length(MLR.panels)
       if strcmp(MLR.panels{iPanel}{1},varargin{1})
	 val = MLR.panels{iPanel}{2};
       end
     end
   end
 case {'callback'}
  % callbackFunctionHandles = viewGet(v,'callback',callbackName)
  % retrieves a cell array of function handles for the callback
  % this is used for parts of the code, like when a base changes,
  % to call external functions that can adjust the behavior of the
  % GUI
  val = {};
  if length(varargin) < 1
    disp(sprintf('(viewGet) Need to specify callbackName'));
    return
  end
  if isfield(MLR,'callbacks')
    % look for the named callback (not case-sensitive)
    callbackNum = find(strcmp(MLR.callbacks{1},lower(varargin{1})));
    if ~isempty(callbackNum)
      % if found, return the callback
      val = MLR.callbacks{2}{callbackNum};
    end
  end

  %% New viewGets go here
  otherwise
    if strcmp(mrGetPref('verbose'),'Yes')
      dispViewGetHelp;
    end
    disp(['Invalid parameter for volume view: ',param]);
end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%   dispViewGetHelp   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function dispViewGetHelp(param)

% open this file and scan the text
fid = fopen(which('viewGet'));
C = textscan(fid,'%s%s','delimiter','{}');
fclose(fid);

collectComments = 0;
commands = {};
for i = 1:length(C{1})
  % collect commens after we have found a case word
  if collectComments & (C{1}{i}(1) == '%')
    if collectComments == 1
      commandComments{length(commands)} = sprintf('   %s',C{1}{i}(3:end));
    else
      commandComments{length(commands)} = sprintf('%s\n   %s',commandComments{length(commands)},C{1}{i}(3:end));
    end
    collectComments = collectComments+1;
  else
    % not a comment, stop collecting comments
    collectComments = 0;
  end
  % get a command
  if (strncmp(C{1}{i},'case',4))
    collectComments = 1;
    commands{end+1} = C{2}{i};
    % get rid of any commands that start with a number (this are spuriuse
    % case clauses from above
    if (length(commands{end})==0) || any(strcmp(commands{end}(1),{'0','1','2','3'}))
      commands = {commands{1:end-1}};
    end
  end
end

% if the user wants help on a particular command
if ~ieNotDefined('param')
  commandNum = [];
  % find the command
  for i = 1:length(commands)
    if ~isempty(findstr(lower(param),commands{i}))
      commandNum(end+1) = i;
    end
  end
  if isempty(commandNum)
    disp(sprintf('(viewGet) No viewGet for %s',param));
  else
    % display all matching commands
    for i = 1:length(commandNum)
      if (commandNum(i) > 0) && (commandNum(i) <= length(commandComments)) && ~isempty(commandComments{commandNum(i)})
        commandString = sprintf('%s %s %s',repmat('=',1,12),commands{commandNum(i)},repmat('=',1,12));
        disp(commandString);
        disp(commandComments{commandNum(i)});
      end
    end
  end
  return
end

% sort everything
[commands sortedIndex] = sort(commands);

maxlen = -inf;
% now print out
for i = 1:length(commands)
  disp(commands{i});
  lens(i) = length(commands{i});
  if length(commandComments)>=sortedIndex(i)
    disp(commandComments{sortedIndex(i)});
  end
end

maxlen = round(median(lens)+4);
nColumns = 6;
disp(sprintf('\n'));
disp('------------------------- All possible parameters ---------------------');
columnNum = 1;
currentChar = '.';
for i = 1:length(commands)
  % make a command with enough spaces
  if ~isempty(commands{i})
    % get this command, remove ',' and replace with /
    thisCommand = fixBadChars(commands{i}(2:end-1),{''',''','/'});
    % if we are switching the start letter, then put a new line
    if thisCommand(1) ~= currentChar
      currentChar = thisCommand(1);
      if columnNum ~= 1
        fprintf(1,sprintf('\n'));
      end
      fprintf(1,'(%s): ',upper(currentChar));
      columnNum = 1;
    elseif columnNum == 1
      fprintf(1,'     ');
    end
    % find out how many columns we need to print out this key word
    thisColumns = ceil((length(thisCommand)+1)/maxlen);
    % now fill up that many columns worth of space
    dispcommand = sprintf(' ');
    dispcommand = repmat(dispcommand,1,thisColumns*maxlen);
    % insert the command name
    dispcommand(1:(length(thisCommand))) = thisCommand;
    % and update what column we are on
    columnNum = columnNum + thisColumns;
    % print the command
    mrDisp(dispcommand);
    % print out new line if we have printed all columns
    if columnNum > nColumns
      fprintf(1,sprintf('\n'));
      columnNum = 1;
    end
  end
end
fprintf(1,sprintf('\n'));
disp('-----------------------------------------------------------------------');

%%%%%%%%%%%%%%%%%%
%%   getGroup   %%
%%%%%%%%%%%%%%%%%%
function g = getGroup(view,varg)

if ieNotDefined('varg')
  g = viewGet(view,'currentGroup');
else
  g = varg{1};
end
if isempty(g)
  g = viewGet(view,'currentGroup');
end
% if group is a string, then convert it to a number
if isstr(g)
  groupName = g;
  g = viewGet(view,'groupNum',g);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%   getScanAndGroup   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function [s g] = getScanAndGroup(view,varg,param,argnum)

if ieNotDefined('argnum'), argnum = 1; end
if ieNotDefined('varg') || (length(varg) < argnum)
  s = viewGet(view,'curScan');
else
  s = varg{argnum};
end
if length(varg) > argnum
  g = varg{argnum+1};
else
  g = viewGet(view,'currentGroup');
end
if isempty(g)
  g = viewGet(view,'currentGroup');
end
if isempty(s)
  s = viewGet(view,'curScan');
end
% if group is a string, then convert it to a number
if isstr(g)
  gnum = viewGet(view,'groupNum',g);
  if isempty(gnum),disp(sprintf('(viewGet:getScanAndGroup) Unknown group: %s', g));end
  g = gnum;
end

%%%%%%%%%%%%%%%%%%%%
%%   getBaseNum   %%
%%%%%%%%%%%%%%%%%%%%
function [b baseVolume] = getBaseNum(view,varg,argnum)

if ieNotDefined('argnum'),argnum = 1;end
if ieNotDefined('varg') || (length(varg) < argnum)
  b = viewGet(view,'currentBase');
else
  b = varg{argnum};
end
if isempty(b)
  b = viewGet(view,'currentBase');
end
if nargout == 2
  baseVolume = viewGet(view,'baseVolume',b);
end

%%%%%%%%%%%%%%%%%%%
%%   getRoiNum   %%
%%%%%%%%%%%%%%%%%%%
function r= getRoiNum(view,varg,argnum)

if ieNotDefined('argnum'), argnum = 1; end
if ieNotDefined('varg') || (length(varg) < argnum)
  r = viewGet(view,'currentROI');
else
  r = varg{argnum};
  if isstr(r)
    r = viewGet(view,'roiNum',r);
    % avoid behaviour where a bad ROI name gives you current ROI instead...
    assert(~isempty(r),'no roi %s found',varg{1});
  end
end
if isempty(r)
  r = viewGet(view,'currentROI');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    viewGetPrependEtc    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fileName = viewGetPrependEtc(view,fileName);

% prepend etc
fileName = fullfile(viewGet(view,'etcDir'),fileName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    viewGetLoadStimFile    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimFile = viewGetLoadStimFile(view,stimFileName)

% prepend etc
stimFileName = fullfile(viewGet(view,'etcDir'),stimFileName);

% load this stimfile
if ~isfile(stimFileName)
  mrErrorDlg(sprintf('(viewGet:viewGetLoadStimfile): Could not find stimfile %s',stimFileName));
else
  stimFile = load(stimFileName);
end
% check to see what type it is, and set the field appropriately
if isfield(stimFile,'mylog')
  stimFile.filetype = 'eventtimes';
  %make sure time and duration vectors are rows
  for iEvent=1:length(stimFile.mylog.stimtimes_s)
    if size(stimFile.mylog.stimtimes_s{iEvent},1)>1
      stimFile.mylog.stimtimes_s{iEvent} = stimFile.mylog.stimtimes_s{iEvent}';
    end
    if isfield(stimFile.mylog,'stimdurations_s') && size(stimFile.mylog.stimdurations_s{iEvent},1)>1
      stimFile.mylog.stimdurations_s{iEvent} = stimFile.mylog.stimdurations_s{iEvent}';
    end
  end
  %move stimNames field
  if isfield(stimFile.mylog,'stimNames')
      stimFile.stimNames = stimFile.mylog.stimNames;
    stimFile.mylog = rmfield(stimFile.mylog,'stimNames');
  end
elseif isfield(stimFile,'stimts')
  stimFile.filetype = 'afni';
elseif isfield(stimFile,'myscreen')
  stimFile.filetype = 'mgl';
elseif isfield(stimFile,'stimvol')
  stimFile.filetype = 'stimvol';
else
  stimFile.filetype = 'unknown';
end
stimFile.filename = stimFileName;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    viewGetPrependPre    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fileName = viewGetPrependPre(view,fileName);

% prepend Pre
fileName = fullfile(viewGet(view,'homeDir'),'Pre',fileName);

%%%%%%%%%%%%%%%%%%%%%%%%
%    viewGetFidInfo    %
%%%%%%%%%%%%%%%%%%%%%%%%
function fidInfo = viewGetFidInfo(view,fileName)

% prepend Pre
fileName = fullfile(viewGet(view,'homeDir'),'Pre',fileName);

% get fidInfo
[xform fidInfo] = fid2xform(fileName);
