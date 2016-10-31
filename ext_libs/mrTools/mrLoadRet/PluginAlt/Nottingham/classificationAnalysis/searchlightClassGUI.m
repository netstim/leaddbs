% surfaceClassGUI.m
%
%      usage: eventRelatedGUI()
%         by: justin gardner
%       date: 04/05/07
%    purpose: 
%
function params = searchlightClassGUI(varargin)

% get the arguments
eval(evalargs(varargin));

if ieNotDefined('defaultParams'),defaultParams = 0;end
if ieNotDefined('thisView'),thisView = newView;end
if ieNotDefined('params'),params = struct;end
if ieNotDefined('groupName'), groupName = viewGet(thisView,'groupName');end;


if ~isfield(params,'groupName') || isempty(params.groupName)
  params.groupName = groupName;
end
if ~isfield(params,'saveName') || isempty(params.saveName)
  params.saveName = 'searchLight';
else
  params.saveName = viewGet(thisView,'analysisName');
end
if ~isfield(params,'scanNum') || isempty(params.scanNum)
  params.scanNum = [];
end
if ~isfield(params,'numberEVs') || isempty(params.numberEvents)
  params.numberEVs = 0;
end
roiMaskMenu = viewGet(thisView,'roiNames');
if isfield(params,'roiMask') && ismember(params.roiMask,roiMaskMenu)
  roiMaskMenu = putOnTopOfList(params.roiMask,roiMaskMenu);
else
end
if ~isfield(params,'radius') || isempty(params.radius)
  params.radius = 1;
end
if ~isfield(params,'sigTest') || isempty(params.sigTest)
  params.sigTest = 1;
end
if ~isfield(params,'fweAdjustment') || isempty(params.fweAdjustment)
  params.fweAdjustment = 0;
end
if ~isfield(params,'fdrAdjustment') || isempty(params.fdrAdjustment)
  params.fdrAdjustment = 0;
end
if ~isfield(params,'advancedMenu') || isempty(params.advancedMenu)
  params.advancedMenu = 0;
end


askForParams = 1;
groupNames = putOnTopOfList(params.groupName,viewGet(thisView,'groupNames'));
while askForParams
% set the parameter string
    paramsInfo = {...
        {'groupName',groupNames,'type=popupmenu','Name of group from which to do searchlight analysis'},...
        {'saveName',params.saveName,'File name to try to save as'},...
        {'roiMask',roiMaskMenu,'ROI to use to mask the analysis, for example an ROI defining the skull-stripped brain, or a GM mask etc'},...
        {'numberEVs',params.numberEVs,'minmax=[0 inf]','incdec=[-1 1]','incdecType=plusMinus','Number of Event Types to be classified. If 0, the number of Event Types will be set to the number of stimulus type in the stimulus file'},...
        {'radius',params.radius,'Radius of searchlight in voxels to use'},...
        {'sigTest',params.sigTest,'type=checkbox','Significance testing of results with Binomial Test'},...
        {'fweAdjustment',params.fweAdjustment,'type=checkbox','contingent=sigTest','Family Wise Error adjustment'},...
        {'fdrAdjustment',params.fdrAdjustment,'type=checkbox','contingent=sigTest','False discovery rate adjustment'},...
%         {'advancedMenu',params.advancedMenu,'type=checkbox','contingent=sigTest','Open menu for advanced statistic options'},...
        };

    % Get parameter values
    if defaultParams
      tempParams = mrParamsDefault(paramsInfo);
    else
      tempParams = mrParamsDialog(paramsInfo,'Searchlight Classification Parameters');
    end

    % if empty user hit cancel
    if isempty(tempParams)
      params = [];
      return;
    end
    params = copyFields(tempParams,params);
    roiMaskMenu = putOnTopOfList(params.roiMask,roiMaskMenu);
    
    while askForParams
        groupNum = viewGet(thisView,'groupnum',params.groupName);
        nScans=viewGet(thisView,'nScans',groupNum);
        if nScans == 1
          params.scanNum = 1;
        elseif ~ieNotDefined('scanList')
          params.scanNum = scanList;
        elseif defaultParams
          params.scanNum = 1:nScans;
        elseif viewGet(thisView,'nScans',groupNum) >1
          scanNums = selectInList(thisView,'scans','Select scans to analyse',params.scanNum,groupNum);
          if isempty(scanNums)
            if size(scanNums,2) %if the top close button has been pressed
              params=[];
              return
            else
              askForParams = 1;
              break;
            end
          else
            params.scanNum = scanNums;
          end
        end
        
        while askForParams
            params.hrfModel='boxCar';params.hrfParams=[];
            [scanParams, params] = getClassVarname(thisView,params,defaultParams);
            if isempty(scanParams)
                if size(scanParams,2) %if the top close button has been pressed
                  params=[];
                  return
                else
                  askForParams = 1;
                  break;
                end
            end
            params.scanParams = scanParams;
            while askForParams
                [scanParams, params] = getClassEVParamsGUI(thisView,params,defaultParams);%getClassEventParamsGUI(thisView,params,defaultParams);
                if isempty(scanParams)
                  if size(scanParams,2) %if the top close button has been pressed
                    params=[];
                    return
                  else
                    askForParams = 1;
                    break;
                  end
                end
                params.scanParams = scanParams;
                askForParams=0;
            end
        end
        if nScans == 1 || ~ieNotDefined('scanList') || defaultParams
          break;
        end
    end
end


% set the scan number
for i = 1:length(params.scanNum)
  params.scanParams{params.scanNum(i)}.scanNum = params.scanNum(i);
end




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % just display parameters
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function retval = dispParams(params)
% 
% paramsInfo = {};
% if isfield(params,'applyFiltering')
%   paramsInfo{end+1} = {'applyFiltering',params.applyFiltering,'editable=0','Whether filtering was applied to the columns of the stimulus convolution matrix'};
% end
% paramsInfo{end+1} = {'analyzedScans',1,'incdec=[-1 1]',sprintf('minmax=[1 %i]',length(params.scanNum)),'editable=0'};
% % get parameters for each scan
% for i = 1:length(params.scanNum)
%   scanNum{i} = params.scanNum(i);
%   description{i} = params.scanParams{scanNum{i}}.description;
%   hdrlen{i} = params.scanParams{scanNum{i}}.hdrlen;
%   preprocess{i} = params.scanParams{scanNum{i}}.preprocess;
%   if isfield(params.scanParams{scanNum{i}},'taskNum')
%     taskNum{i} = params.scanParams{scanNum{i}}.taskNum;
%   end
%   if isfield(params.scanParams{scanNum{i}},'phaseNum')
%     phaseNum{i} = params.scanParams{scanNum{i}}.phaseNum;
%   end
%   if isfield(params.scanParams{scanNum{i}},'segmentNum')
%     segmentNum{i} = params.scanParams{scanNum{i}}.segmentNum;
%   end
%   if isfield(params.scanParams{scanNum{i}},'varname')
%     varname{i} = params.scanParams{scanNum{i}}.varname;
%   end
% end
% paramsInfo{end+1} = {'scanNum',scanNum,'group=analyzedScans','type=numeric','editable=0','scan number'};
% paramsInfo{end+1} = {'tseriesFile',params.tseriesFile,'group=analyzedScans','type=string','editable=0','Name of timeseries that was analyzed'};
% paramsInfo{end+1} = {'description',description,'group=analyzedScans','type=string','editable=0','Description of the analysis'};
% paramsInfo{end+1} = {'eventLength',eventLength,'group=analyzedScans','type=numeric','editable=0','Length of response in seconds to calculate'};
% paramsInfo{end+1} = {'preprocess',preprocess,'group=analyzedScans','type=string','editable=0','String of extra commands for preprocessing. Normally you will not need to set anything here, but this allows you to do corrections to the stimvols that are calculated so that you can modify the analysis. (see wiki for details)'};
% if ~ieNotDefined('taskNum')
%   paramsInfo{end+1} = {'taskNum',taskNum,'group=analyzedScans','type=numeric','editable=0','The task you want to use'};
% end
% if ~ieNotDefined('phaseNum')
%   paramsInfo{end+1} = {'phaseNum',phaseNum,'group=analyzedScans','type=numeric','editable=0','The phase of the task you want to use'};
% end
% if ~ieNotDefined('segmentNum')
%   paramsInfo{end+1} = {'segmentNum',segmentNum,'group=analyzedScans','type=numeric','editable=0','The segment of the task you want to use'};
% end
% if ~ieNotDefined('varname')
%   paramsInfo{end+1} = {'varname',varname,'group=analyzedScans','type=string','editable=0','The variable that was analyzed'};
% end
% retval = mrParamsDialog(paramsInfo,'Event Related parameters');
