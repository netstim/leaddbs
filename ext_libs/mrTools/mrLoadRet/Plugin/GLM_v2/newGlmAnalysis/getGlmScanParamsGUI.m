% getScanParamsGUI.m
%
%        $Id$
%      usage: params = getScanParamsGUI(thisView,params,useDefault)
%         by: julien besle, extracted from eventRelatedGLMGUI so taht it can be called from eventRelatedGUI eventually
%       date: 23/03/2010
%    purpose: return scan specific parameters for GLM analysis
%

function [scanParams, params] = getGlmScanParamsGUI(thisView,params,useDefault)

keepAsking = 1;
groupNum = viewGet(thisView,'groupNum',params.groupName);
nScans = viewGet(thisView,'nScans',groupNum);
if ~isfield(params,'scanNum') || isempty(params.scanNum)
  params.scanNum = 1:nScans;
end
if isfield(params,'scanParams') 
  if length(params.scanParams)>=nScans
    scanParams = params.scanParams(1:nScans);
  else
    scanParams = cell(1,nScans);
    scanParams(1:length(params.scanParams))=params.scanParams;
  end
else
   % make the output as long as the number of scans
   scanParams = cell(1,nScans);
end

while keepAsking
  % check for stimfile, and if it is mgl/type then ask the
  % user which variable they want to do the anlysis on
  for iScan = params.scanNum
   % get scan info and description
    if ~isfield(scanParams{iScan},'scan')
       scanParams{iScan}.scan = sprintf('%i: %s',iScan,viewGet(thisView,'description',iScan,groupNum));
    end
    if ~isfield(scanParams{iScan},'description')
       scanParams{iScan}.description = sprintf('GLM analysis of %s: %i',params.groupName,iScan);
    end
    if ~isfield(scanParams{iScan},'preprocess')
       scanParams{iScan}.preprocess = '';
    end
    if ~isfield(scanParams{iScan},'subsetBox') || ~strcmp(params.analysisVolume,'Subset box')
       scanDims = viewGet(thisView,'dims',iScan,groupNum);
       scanParams{iScan}.subsetBox = ['[1 ' num2str(scanDims(1)) ';1 ' num2str(scanDims(2)) ';1 ' num2str(scanDims(3)) ']'];
    end
    if ~isfield(scanParams{iScan},'sameForNextScans') || isempty(scanParams{iScan}.sameForNextScans)
       scanParams{iScan}.sameForNextScans = iScan == params.scanNum(1);
    end
    if fieldIsNotDefined(scanParams{iScan},'highpassDesign')
       scanParams{iScan}.highpassDesign = 1;
    end
    if fieldIsNotDefined(scanParams{iScan},'projectOutGlobalComponent')
       scanParams{iScan}.projectOutGlobalComponent = 1;
    end

 % Standard parameters to set
    
    subsetBoxVisibleOption = 'visible=1';
    subsetBoxEditableOption = 'editable=1';
    if ismember(params.analysisVolume,{'Loaded ROI(s)' 'Visible ROI(s)'})
      subsetBoxVisibleOption = 'visible=0';
    elseif strcmp(params.analysisVolume,'Whole volume')
      subsetBoxEditableOption = 'editable=0';
    end
    subsetBoxHelp = 'Subset of voxels on which to perform the analysis, in the form [X1 X2;Y1 Y2;Z1 Z2] (Zs are optional).'; 
    if params.covCorrection
      subsetBoxHelp = sprintf('%s In order to estimate the noise covariance matrix for voxels of interest, this box should also include margins of %d voxels of no interest on both sides of the X and Y dimensions.',subsetBoxHelp,floor(params.covEstimationAreaSize/2));
    end
    
    d= loadScan(thisView, iScan, groupNum, 0);
    
    paramsInfo = {...
    {'scan',scanParams{iScan}.scan,'type=statictext','Description of scan to set parameters for (not editable)'},...
    {'description',scanParams{iScan}.description,'Event related analysis of [x...x]','Description of the analysis'}...
    {'preprocess',scanParams{iScan}.preprocess,'String of extra commands for preprocessing. Normally you will not need to set anything here, but this allows you to do corrections to the stimvols that are calculated so that you can modify the analysis. (see wiki for details)',...
                  'callback',{@tryPreProcess,d},'passParams=1'}...
    {'subsetBox', scanParams{iScan}.subsetBox, subsetBoxVisibleOption, subsetBoxEditableOption, subsetBoxHelp};
         };

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
        task = cellArray(stimfile{1}.task,2);
        taskNumList = num2cell(1:length(task));
        if isfield(scanParams{iScan},'taskNum')
          taskNumList = putOnTopOfList(scanParams{iScan}.taskNum,taskNumList);
        end
        % if there is more than one task, then ask the user for that
        if length(task)>1
          paramsInfo{end+1} = {'taskNum',taskNumList,'The task you want to use'};
        end
        % if there are multiple phases, then ask for that
        maxPhaseNum = 0;
        maxSegNum = 0;
        phaseNumList = {};
        for tnum = 1:length(task)
          phaseNumList{tnum} = num2cell(1:length(task{tnum}));
          maxPhaseNum = max(maxPhaseNum,length(task{tnum}));
          % if there are multiple _segments_, then ask for that
          for pnum = 1:length(task{tnum})
            segNum{tnum}{pnum} = num2cell(1:length(task{tnum}{pnum}.segmin));
            maxSegNum = max(maxSegNum,length(segNum{tnum}{pnum}));
          end
        end
        if maxPhaseNum > 1
          if ~isfield(scanParams{iScan},'phaseNum')
            scanParams{iScan}.phaseNum = phaseNumList{taskNumList{1}}{1};
          end
          phaseNumList{taskNumList{1}} = putOnTopOfList(scanParams{iScan}.phaseNum, phaseNumList{taskNumList{1}});
          if length(task) == 1
            paramsInfo{end+1} = {'phaseNum',phaseNumList{1},'The phase of the task you want to use'};
          else
            paramsInfo{end+1} = {'phaseNum',phaseNumList,'The phase of the task you want to use','contingent=taskNum'};
          end
        end

         % if there is more than one segment in any of the phases, ask the user to specify
         % should add some error checking.
        if maxSegNum > 1
          if ~isfield(scanParams{iScan},'segmentNum')
             scanParams{iScan}.segmentNum = 1;
          end
          paramsInfo{end+1} = {'segmentNum',scanParams{iScan}.segmentNum,'The segment of the task you want to use','incdec=[-1 1]','incdecType=plusMinus'};
        end

        %set up to get the variable name from the user
        if ~isfield(scanParams{iScan},'varname')
          scanParams{iScan}.varname = varnames{1};
        end
        if iscell(scanParams{iScan}.varname)
            %convert cell array to string
            varname = '{';
            for iCell = 1:length(scanParams{iScan}.varname)
                if iscell(scanParams{iScan}.varname{iCell})
                    varname = [varname '{'];
                    for jCell = 1:length(scanParams{iScan}.varname{iCell})
                        varname = [varname '''' scanParams{iScan}.varname{iCell}{jCell} ''',',];
                    end
                    varname = [varname(1:end-1) '},'];
                else
                    varname = [varname '''' scanParams{iScan}.varname{iCell} ''',',];
                end
            end
            varname = [varname(1:end-1) '}'];
        else
            varname = scanParams{iScan}.varname;
        end
        paramsInfo{end+1} ={'varname',varname,sprintf('Analysis variables: %s',varnamesStr)};
      end
    end
    
    concatenationOptions = 0;
    if isfield(d,'concatInfo') %in case of concatenated files
      % ask if hipass filter should be applied to the design matrix
      if isfield(d.concatInfo,'hipassfilter') && ~isempty(d.concatInfo.hipassfilter{1})
        paramsInfo = [paramsInfo {...
          {'highpassDesign',scanParams{iScan}.highpassDesign, 'type=checkbox','For concatenated scans, applies to the design matrix the high-pass filter that has been applied to the individual scans.'},...
           }];
        concatenationOptions = 1;
      end
      % ask if  mean vector should be projected out the
      if isfield(d.concatInfo,'projection') && ~isempty(d.concatInfo.projection{1})
        paramsInfo = [paramsInfo {...
          {'projectOutGlobalComponent',scanParams{iScan}.highpassDesign, 'type=checkbox','For concatenated scans, projects out the global activity estimated on an ROI from the design matrix.'},...
           }];
        concatenationOptions = 1;
      end
    end
     
     % give the option to use the same variable for remaining scans
    if (iScan ~= params.scanNum(end)) && (length(params.scanNum)>1)
      paramsInfo{end+1} = {'sameForNextScans',scanParams{iScan}.sameForNextScans,'type=checkbox','Use the same parameters for all scans'};
    end

    %%%%%%%%%%%%%%%%%%%%%%%
    % now we have all the dialog information, ask the user to set parameters
    if useDefault || (isempty(strfind(stimfile{1}.filetype,'mgl')) && ~strcmp(params.analysisVolume,'Subset box') && ~concatenationOptions)
       tempParams = mrParamsDefault(paramsInfo);
    else
       tempParams = mrParamsDialog(paramsInfo,'Set Scan Parameters');
    end

    % user hit cancel
    if isempty(tempParams)
       scanParams = tempParams;
       return
    end
    
    
    subsetBox = eval(tempParams.subsetBox);
    
    %various controls and variables settings
    if params.covCorrection && strcmp(params.analysisVolume,'Subset box') any(diff(subsetBox(1:2,:),1,2)<=params.covEstimationAreaSize)
      mrWarnDlg('(getScanParamsGUI) subset box size is too small for covariance estimation area size','Yes');
      scanParams{iScan} = [];
    else
      keepAsking = 0;

      % check if the varname is a cell array, then convert to a cell array
      % instead of a string this is so that the user can specify a variable
      % name like {{'varname'}}
      if (isfield(tempParams,'varname') &&...
         ischar(tempParams.varname) && ...
         (length(tempParams.varname) > 1) && ...
         (tempParams.varname(1) == '{'))
         tempParams.varname = eval(tempParams.varname);
      end

      scanParams{iScan} = mrParamsCopyFields(tempParams,scanParams{iScan});
      
      % if sameForNextScans is set, copy all parameters into remaining scans and break out of loop
      if isfield(tempParams,'sameForNextScans') && tempParams.sameForNextScans
         for jScan = params.scanNum(find(params.scanNum>iScan,1,'first'):end)
            % set the other scans params to the same as this one
            scanParams{jScan} = mrParamsCopyFields(tempParams,scanParams{jScan});
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
      
    end
    
    if keepAsking && useDefault %there were incompatible parameters but this is the script mode (no GUI)
      scanParams = [];
      return;
    end
  end
  
  %if analysis on whole volume, need to set subset box to scan dimensions 
  if ~strcmp(params.analysisVolume,'Subset box')
    for iScan = params.scanNum(1:end)
      scanDims = viewGet(thisView,'dims',iScan,groupNum);
      scanParams{iScan}.subsetBox = ['[1 ' num2str(scanDims(1)) ';1 ' num2str(scanDims(2)) ';1 ' num2str(scanDims(3)) ']'];
    end
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% tryPreProcess %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function preprocess = tryPreProcess(scanParams,d)

% Try the pre-processing function
preProcessFailure = 0;
if ~isempty(scanParams.preprocess)
  d = getStimvol(d,scanParams);
  [d, preProcessFailure] = eventRelatedPreProcess(d,scanParams.preprocess);
end

if preProcessFailure
  mrWarnDlg(['(getScanParamsGUI) There was a problem running pre-processing function ' scanParams.preprocess],'Yes');
  preprocess = [];
end

