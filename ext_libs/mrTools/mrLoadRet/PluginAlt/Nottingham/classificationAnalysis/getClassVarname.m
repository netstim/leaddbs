% getEventRelatedVarname.m
%
%      usage: getEventRelatedVarname(view,<groupNum>,<scanNums>,<useDefault>)
%         by: justin gardner
%       date: 06/25/08
%    purpose: function to get the variable name that the user wants
%             to do the event related analysis on, puts up a gui
%
function [scanParams, params] = getClassVarname(thisView,params,useDefault)

keepAsking = 1;
groupNum = viewGet(thisView,'groupNum',params.groupName);
nScans = viewGet(thisView,'nScans',groupNum);
if ~isfield(params,'scanNum') || isempty(params.scanNum)
  params.scanNum = 1:nScans;
end
if isfield(params,'scanParams') && length(params.scanParams)==nScans
   scanParams = params.scanParams;
else
   % make the output as long as the number of scans
   scanParams = cell(1,nScans);
end

while keepAsking
    for iScan=params.scanNum
        tr = viewGet(thisView,'framePeriod',iScan,groupNum);
        if ~isfield(scanParams{iScan},'scan')
           scanParams{iScan}.scan = sprintf('%i: %s',iScan,viewGet(thisView,'description',iScan,groupNum));
        end
        if ~isfield(scanParams{iScan},'description')
           scanParams{iScan}.description = sprintf('ROI-classification analysis of %s: %i',params.groupName,iScan);
        end
        if ~isfield(scanParams{iScan},'preprocess')
           scanParams{iScan}.preprocess = '';
        end
        if ~isfield(scanParams{iScan},'sameForNextScans') || isempty(scanParams{iScan}.sameForNextScans)
            scanParams{iScan}.sameForNextScans = iScan == params.scanNum(1);
        end
        paramsInfo = {...
            {'scan',scanParams{iScan}.scan,'type=statictext','Description of scan to set parameters for (not editable)'},...
            {'description',scanParams{iScan}.description,'ROI-classification analysis of [x...x]','Description of the analysis'}...
            {'preprocess',scanParams{iScan}.preprocess,'String of extra commands for preprocessing. Normally you will not need to set anything here, but this allows you to do corrections to the stimvols that are calculated so that you can modify the analysis. (see wiki for details)',...
                  'callback',{@tryPreProcess,loadScan(thisView, iScan, groupNum, 0)},'passParams=1'}...
                  };
              
              % Timing parameters
    % make sure we are running on a set with a stimfile
    stimfile = viewGet(thisView,'stimfile',iScan,groupNum);
    
    if isempty(stimfile)
     mrMsgBox(sprintf('No associated stimfile with scan %i in group %s',iScan,params.groupName));
     scanParams = [];
     return
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % see if we have a stimfile from mgl, in which case we should
    % ask the user what the variable name is that they want to use for the analysis
    if strfind(stimfile{1}.filetype,'mgl')

      % check to see what style this is, if the task variable does
      % not have a segmentTrace then it mus be an old style, in which
      % we used channels
      task = cellArray(stimfile{1}.task,2);
      if isfield(stimfile{1}.myscreen,'traces') && ~isfield(task{1}{1},'segmentTrace')
        % this is the old style, get the stimtrace number
        paramsInfo{end+1} = {'stimtrace',stimfile{1}.myscreen.stimtrace,'the trace number that contains the stimulus','incdec=[-1 1]','incdecType=plusMinus',sprintf('minmax=[%i %i]',stimfile{1}.myscreen.stimtrace,size(stimfile{1}.myscreen.traces,1))};
      else
        if exist('getTaskVarnames') ~= 2
          mrErrorDlg('(eventRelatedGUI) MGL function getTaskVarnames is not in path. You must have mgl in the path to extract stimulus timing information from an mgl stim file');
        end
        % this is the new tyle, ask for a variable name
        [varnames varnamesStr] = getTaskVarnames(stimfile{1}.task);
        % if there is more than one task, then ask the user for that
        task = cellArray(stimfile{1}.task,2);
        if length(task)>1
          paramsInfo{end+1} = {'taskNum',num2cell(1:length(task)),'The task you want to use'};
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
            paramsInfo{end+1} = {'phaseNum',phaseNum{1},'The phase of the task you want to use'};
          else
            paramsInfo{end+1} = {'phaseNum',phaseNum,'The phase of the task you want to use','contingent=taskNum'};
          end
        end

         % if there is more than one segment in any of the phases, ask the user to specify
         % should add some error checking.
        if maxSegNum > 1
          paramsInfo{end+1} = {'segmentNum',1,'The segment of the task you want to use','incdec=[-1 1]','incdecType=plusMinus'};
        end

         %set up to get the variable name from the user
%         paramsInfo{end+1} ={'varname',varnames{1},sprintf('Analysis variables: %s',varnamesStr),'type=popupmenu'};
        paramsInfo{end+1} ={'varname',varnames,sprintf('Analysis variables: %s',varnamesStr),'type=popupmenu'};
        if ~isfield(scanParams{iScan},'stimDuration') || strcmp(params.hrfModel,'hrfDeconvolution')
            scanParams{iScan}.stimDuration = num2str(tr);
        elseif isnumeric(scanParams{iScan}.stimDuration)
            scanParams{iScan}.stimDuration = num2str(scanParams{iScan}.stimDuration);
        end
%          paramsInfo{end+1} ={'stimDuration', scanParams{iScan}.stimDuration, 'duration of stimulation/event (seconds, min=0.01s), a boxcar function that is convolved with hrf.'};
      end
    elseif strfind(stimfile{1}.filetype,'eventtimes')  && ~isfield(scanParams{iScan},'stimDuration') && isfield(stimfile{1}.mylog,'stimdurations_s')
          scanParams{iScan}.stimDuration = 'fromFile';
    end
    
    if ~isfield(scanParams{iScan},'stimDuration')
       scanParams{iScan}.stimDuration = num2str(tr);
    elseif isnumeric(scanParams{iScan}.stimDuration)
      scanParams{iScan}.stimDuration = num2str(scanParams{iScan}.stimDuration);
    end
    paramsInfo{end+1} ={'stimDuration', scanParams{iScan}.stimDuration, 'duration of stimulation/event (seconds, min=0.01s), a boxcar function that is convolved with hrf.'};
    
     % give the option to use the same variable for remaining scans
        if (iScan ~= params.scanNum(end)) && (length(params.scanNum)>1)
            paramsInfo{end+1} = {'sameForNextScans',scanParams{iScan}.sameForNextScans,'type=checkbox','Use the same parameters for all scans'};
        end
        
        % now we have all the dialog information, ask the user to set parameters
    if useDefault
       tempParams = mrParamsDefault(paramsInfo);
    else
       tempParams = mrParamsDialog(paramsInfo,'Set Scan Parameters');
    end

    % user hit cancel
    if isempty(tempParams)
       scanParams = tempParams;
       return
    end
    tempParams.stimDuration= str2num(tempParams.stimDuration)
    
    scanParams{iScan} = mrParamsCopyFields(tempParams,scanParams{iScan});
    
    keepAsking = 0;

      % check if the varname is a cell array, then convert to a cell array
      % instead of a string this is so that the user can specify a variable
      % name like {{'varname'}}
      if (isfield(scanParams{iScan},'varname') &&...
         ischar(scanParams{iScan}.varname) && ...
         (length(scanParams{iScan}.varname) > 1) && ...
         (scanParams{iScan}.varname(1) == '{'))
         scanParams{iScan}.varname = eval(scanParams{iScan}.varname);
      end
%       scanParams{iScan}.stimDuration = scanParams{iScan}.eventLength*2;
      % if sameForNextScans is set, copy all parameters into remaining scans and break out of loop
      if isfield(scanParams{iScan},'sameForNextScans') && ...
         scanParams{iScan}.sameForNextScans
         for jScan = params.scanNum(find(params.scanNum>iScan,1,'first'):end)
            % set the other scans params to the same as this one
            scanParams{jScan} = scanParams{iScan};
            % change the description and scan fields appropriately for this scan num
            scanParams{jScan}.scan = viewGet(thisView,'description',jScan,groupNum);
            description = scanParams{iScan}.description;
            groupNameLoc = strfind(description,params.groupName);
            if ~isempty(groupNameLoc)
               description = sprintf('%s%s: %i',description(1:groupNameLoc(1)),params.groupName,jScan);
            end
            scanParams{jScan}.description = description;
         end
         break;
      %else if the next scan params are empty, copy those from this scan  
      elseif (iScan ~= params.scanNum(end)) 
        nextScan = params.scanNum(find(params.scanNum>iScan,1,'first'));
        if isempty(scanParams{nextScan})
          scanParams{nextScan} = scanParams{iScan};
          scanParams{nextScan} = mrParamsRemoveField(scanParams{nextScan},'scan');
          scanParams{nextScan} = mrParamsRemoveField(scanParams{nextScan},'description');
        end
      end
    if keepAsking && useDefault %there were incompatible parameters but this is the script mode (no GUI)
      scanParams = [];
      return;
    end
    
    end
end


return


function preprocess = tryPreProcess(scanParams,d)

% Try the pre-processing function
preProcessFailure = 0;
if ~isempty(scanParams.preprocess)
  [d, preProcessFailure] = eventRelatedPreProcess(d,scanParams.preprocess);
end

if preProcessFailure
  mrWarnDlg(['(getScanParamsGUI) There was a problem running pre-processing function ' scanParams.preprocess],'Yes');
  preprocess = [];
end

return





%%


% get default variables if necessary
if ieNotDefined('useDefault'),useDefault = 0;end
if ieNotDefined('groupNum'),groupNum = viewGet(view,'curGroup');end
groupName = viewGet(view,'groupName',groupNum);
if ieNotDefined('scanNums'),scanNums = viewGet(view,'curScan');end
if ieNotDefined('eventLength'),eventLength = 1;end
% if ieNotDefined('radiusSize'),radiusSize = 1;end
% make the output as long as the number of scans
scanParams = cell(1,viewGet(view,'nScans',groupNum));

% check for stimfile, and if it is mgl/type then ask the
% user which variable they want to do the anlysis on
for scanNum = 1:length(scanNums)
  % get scan and default description
  scanInfo = sprintf('%i: %s',scanNums(scanNum),viewGet(view,'description',scanNums(scanNum)));
  description = sprintf('Searchlight Classification analysis of %s: %i',groupName,scanNums(scanNum));
  % standard parameters to set
  taskVarParams = {...
      {'scan',scanInfo,'type=statictext','Description of scan to set parameters for (not editable)'},...
      {'description',description,'Event related analysis of [x...x]','Description of the analysis'}...
      {'eventLength',eventLength,'Length of your stimulus in TRs'}...
      {'hdLag',0,'Estimated hemodynamic lag in TRs'}...
      {'averageEvent',1,'type=checkbox','Whether to average across the TRs for each stimulus'}...
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
      taskVarParams{end+1} ={'varname',varnames,sprintf('Analysis variables: %s',varnamesStr),'type=popupmenu'};
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

