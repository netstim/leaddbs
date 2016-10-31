% getEventRelatedVarname.m
%
%      usage: getEventRelatedVarname(view,<groupNum>,<scanNums>,<useDefault>)
%         by: justin gardner
%       date: 06/25/08
%    purpose: function to get the variable name that the user wants
%             to do the event related analysis on, puts up a gui
%
function scanParams = getEventRelatedVarname(view,groupNum,scanNums,useDefault,varargin)

% check arguments
if ~any(nargin == [1 2 3 4 5])
  help getEventRelatedVarname
  return
end

eval(evalargs(varargin,[],[],{'hdrlen'}));

% get default variables if necessary
if ieNotDefined('useDefault'),useDefault = 0;end
if ieNotDefined('groupNum'),groupNum = viewGet(view,'curGroup');end
groupName = viewGet(view,'groupName',groupNum);
if ieNotDefined('scanNums'),scanNums = viewGet(view,'curScan');end
if ieNotDefined('hdrlen'),hdrlen = 25;end
% make the output as long as the number of scans
scanParams = cell(1,viewGet(view,'nScans',groupNum));

% check for stimfile, and if it is mgl/type then ask the
% user which variable they want to do the anlysis on
for scanNum = 1:length(scanNums)
  % get scan and default description
  scanInfo = sprintf('%i: %s',scanNums(scanNum),viewGet(view,'description',scanNums(scanNum)));
  description = sprintf('Event related analysis of %s: %i',groupName,scanNums(scanNum));
  % standard parameters to set
  taskVarParams = {...
      {'scan',scanInfo,'type=statictext','Description of scan to set parameters for (not editable)'},...
      {'description',description,'Event related analysis of [x...x]','Description of the analysis'}...
      {'hdrlen',hdrlen,'Length of response in seconds to calculate'}...
      {'preprocess','','String of extra commands for preprocessing. Normally you will not need to set anything here, but this allows you to do corrections to the stimvols that are calculated so that you can modify the analysis. (see wiki for details)'}...
		  };

  % make sure we are running on a set with a stimfile
  stimfile = viewGet(view,'stimfile',scanNums(scanNum));
  
  if isempty(stimfile)
    mrMsgBox(sprintf('No associated stimfile with scan %i in group %s',scanNums(scanNum),groupName));
    scanParams = [];
    return
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % see if we have a stimfile from mgl, in which case we should
  % ask the user what the variable name is that they want ot use for the analysis
  if strfind(stimfile{1}.filetype,'mgl')

    % check to see what style this is, if the task variable does
    % not have a segmentTrace then it mus be an old style, in which
    % we used channels
    task = cellArray(stimfile{1}.task,2);
    if isfield(stimfile{1}.myscreen,'traces') && ~isfield(task{1}{1},'segmentTrace')
      % this is the old style, get the stimtrace number
      taskVarParams{end+1} = {'stimtrace',stimfile{1}.myscreen.stimtrace,'the trace number that contains the stimulus','incdec=[-1 1]',sprintf('minmax=[%i %i]',stimfile{1}.myscreen.stimtrace,size(stimfile{1}.myscreen.traces,1))};
    else
      if exist('getTaskVarnames') ~= 2
	mrErrorDlg('(eventRelatedGUI) MGL function getTaskVarnames is not in path. You must have mgl in the path to extract stimulus timing information from an mgl stim file');
      end
      % this is the new tyle, ask for a variable name
      [varnames varnamesStr] = getTaskVarnames(stimfile{1}.task);
      if isempty(varnames)
	mrErrorDlg('(eventRelatedGUI) No varnames found in stimfile');
      end
      % if there is more than one task, then ask the user for that
      task = cellArray(stimfile{1}.task,2);
      if length(task)>1
	taskVarParams{end+1} = {'taskNum',num2cell(1:length(task)),'The task you want to use'};
      end
      % if there are multiple phases, then ask for that
      maxPhaseNum = 0;
      maxSegNum = 0;
      for tnum = 1:length(task)
	phaseNum{tnum} = num2cell(1:length(task{tnum}));
	maxPhaseNum = max(maxPhaseNum,length(task{tnum}));
	% if there are multiple _segments_, then ask for that
	for pnum = 1:length(task{tnum})
	  segNum{tnum}{pnum} = num2cell(1:length(task{tnum}{pnum}.segmin));
	  maxSegNum = max(maxSegNum,length(segNum{tnum}{pnum}));
	end
      end
      if maxPhaseNum > 1
	if length(task) == 1
	  taskVarParams{end+1} = {'phaseNum',phaseNum{1},'The phase of the task you want to use'};
	else
	  taskVarParams{end+1} = {'phaseNum',phaseNum,'The phase of the task you want to use','contingent=taskNum'};
	end
      end
      
      % if there is more than one segement in any of the phases, ask the user to specify
      % should add some error checking.
      if maxSegNum > 1
	  taskVarParams{end+1} = {'segmentNum',1,'The segment of the task you want to use','incdec=[-1 1]'};
      end
      
      % set up to get the variable name from the user
      taskVarParams{end+1} ={'varname',varnames{1},sprintf('Analysis variables: %s',varnamesStr)};
    end
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  % give the option to use the same variable for all
  if (scanNum == 1) && (length(scanNums)>1)
    taskVarParams{end+1} = {'sameForAll',1,'type=checkbox','Use the same variable name for all analyses'};
  end
  %%%%%%%%%%%%%%%%%%%%%%%
  % now we have all the dialog information, ask the user to set parameters
  if useDefault
    scanParams{scanNums(scanNum)} = mrParamsDefault(taskVarParams);
  else
    scanParams{scanNums(scanNum)} = mrParamsDialog(taskVarParams);
  end
  % user hit cancel
  if isempty(scanParams{scanNums(scanNum)})
    scanParams = [];
    return
  end
  %%%%%%%%%%%%%%%%%%%%%%%
    
  % check if the varname is a cell array, then convert to a cell array
  % instead of a string this is so that the user can specify a variable
  % name like {{'varname'}}
  if (isfield(scanParams{scanNums(scanNum)},'varname') &&...
      ischar(scanParams{scanNums(scanNum)}.varname) && ...
      (length(scanParams{scanNums(scanNum)}.varname) > 1) && ...
      (scanParams{scanNums(scanNum)}.varname(1) == '{'))
    scanParams{scanNums(scanNum)}.varname = eval(scanParams{scanNums(scanNum)}.varname);
  end

  % if sameForAll is set, copy all parameters into all scans and break out of loop
  if isfield(scanParams{scanNums(scanNum)},'sameForAll') && ...
	scanParams{scanNums(scanNum)}.sameForAll
    for i = 2:length(scanNums)
      % set the other scans params to the same as this one
      scanParams{scanNums(i)} = scanParams{scanNums(1)};
      % change the description field appropriately for this scan num
      description = scanParams{scanNums(1)}.description;
      groupNameLoc = strfind(description,groupName);
      if ~isempty(groupNameLoc)
	description = sprintf('%s%s: %i',description(1:groupNameLoc(1)),groupName,scanNums(i));
      end
      scanParams{scanNums(i)}.description = description;
    end
    break
  end
  taskVarParams = {};
end

