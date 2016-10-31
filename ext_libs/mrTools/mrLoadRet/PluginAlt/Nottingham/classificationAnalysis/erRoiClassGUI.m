% eventRelatedGUI.m
%
%      usage: eventRelatedGUI()
%         by: justin gardner
%       date: 04/05/07
%    purpose: 
%
function params = erRoiClassGUI(varargin)

% get the arguments
eval(evalargs(varargin));

if ieNotDefined('defaultParams'),defaultParams = 0;end
if ieNotDefined('thisView'),thisView = newView;end
if ieNotDefined('params'),params = struct;end
if ieNotDefined('groupName'), groupName = viewGet(thisView,'groupName');end;

if ~isfield(params,'hrfParams')
  params.hrfParams = struct;
end
if ~isfield(params,'groupName') || isempty(params.groupName)
  params.groupName = groupName;
end
if ~isfield(params,'saveName') || isempty(params.saveName)
  params.saveName = 'roiClass';
end
if ~isfield(params,'scanNum') || isempty(params.scanNum)
  params.scanNum = [];
end
if ~isfield(params,'numberEVs') || isempty(params.numberEVs)
  params.numberEVs = 0;
end
if ~isfield(params,'diagLinear') || isempty(params.diagLinear)
  params.diagLinear = 1;
end
if ~isfield(params,'Linear') || isempty(params.Linear)
  params.Linear = 0;
end
if ~isfield(params,'SVM') || isempty(params.SVM)
  params.SVM = 0;
end
if ~isfield(params,'sigTest') || isempty(params.sigTest)
  params.sigTest = 0;
end
if ~isfield(params,'nonParaTest') || isempty(params.nonParaTest)
  params.nonParaTest = 0;
end
if ~isfield(params,'numShuff') || isempty(params.numShuff)
  params.numShuff = 0;
end
if ~isfield(params,'selectVox') || isempty(params.selectVox)
  params.selectVox = 0;
end
if ~isfield(params,'hrfModel') || isempty(params.hrfModel)
  params.hrfModel = 'boxCar';
end
if ~isfield(params,'patternMethod') || isempty(params.patternMethod)
  params.patternMethod = 'onePerRun';
end
% if ~isfield(params,'advancedMenu') || isempty(params.advancedMenu)
%   params.advancedMenu = 0;
% end


askForParams = 1;
groupNames = putOnTopOfList(params.groupName,viewGet(thisView,'groupNames'));
hrfModelMenu = putOnTopOfList(params.hrfModel,{'boxCar','hrfDoubleGamma'});
patternMethodMenu = putOnTopOfList(params.patternMethod,{'onePerRun','onePerTrial'});
while askForParams
    % set the parameter string
    paramsInfo = {...
        {'groupName',groupNames,'type=popupmenu','Name of group from which to do eventRelated analysis'},...
        {'saveName',params.saveName,'File name to try to save as'},...
        {'numberEVs',params.numberEVs,'minmax=[0 inf]','incdec=[-1 1]','incdecType=plusMinus','Number of Event Types to be classified. If 0, the number of Event Types will be set to the number of stimulus type in the stimulus file'},...
        {'patternMethod',patternMethodMenu,'type=popupmenu','Method for calculating patterns, either one per run (i.e. Kamitani & Tong) or one per trial (i.e. Kourtzi et al)'},...
        {'hrfModel',hrfModelMenu,'type=popupmenu','Name of the function that defines the Hemodynamic Response Function model that will be convolved with the design matrix. Default = averaging boxCar',},...
        {'diagLinear',params.diagLinear,'type=checkbox','Use diagLinear classification'}...
        {'Linear',params.Linear,'type=checkbox','Use Linear classification'}...
        {'SVM',params.SVM,'type=checkbox','Use svm classification'}...
        {'sigTest',params.sigTest,'type=checkbox','Significance Testing with a Binomial'}...
        {'nonParaTest',params.nonParaTest,'type=checkbox','Significance Testing using non-paramatric shuffling'}...
        {'numShuff',params.numShuff,'contingent=nonParaTest','Number of label reshuffles to do'}...
        {'selectVox',params.selectVox,'How many voxels to use in each ROI, will require selection of an statistical overlay to sort by. Leave 0 to use all voxels in each ROI'},...
    %     {'inplaceConcat',0,'type=checkbox','Concatenate all data and design matrices in memory. This runs a differrent processing stream (ask Farshad for details). If you are using a Concatenation time series do not check this.'},...
    %     {'applyFiltering',1,'type=checkbox','Apply the same filtering that was used before concatenation to the design matrix. This will apply the hipass filter as well as the projection analsis if they have been done by concatTSeries'},...
    };

    % Get parameter values
    if defaultParams
      tempParams = mrParamsDefault(paramsInfo);
    else
      tempParams = mrParamsDialog(paramsInfo,'ROI Classification Parameters');
    end

    % if empty user hit cancel
    if isempty(tempParams)
      params = [];
      return;
    end
    if ~any([tempParams.diagLinear tempParams.Linear tempParams.SVM])
        mrWarnDlg('No classification method selected, using diagLinear by default');
       tempParams.diagLinear = 1;
    end
    hrfModelMenu = putOnTopOfList(tempParams.hrfModel,hrfModelMenu);

    params = copyFields(tempParams,params);
    while askForParams       % get hrf model parameters
      %here we assume that all scans in this group have the same framePeriod
      hrfParams = feval(params.hrfModel,params.hrfParams,viewGet(thisView,'framePeriod',1,viewGet(thisView,'groupNum',params.groupName)),0,defaultParams);
      % if empty user hit cancel, go back
      if isempty(hrfParams)
        if size(hrfParams,2) %if the top close button has been pressed
          params=[];
          return
        else
          break;
        end
      end
      params.hrfParams = hrfParams;
      if ~isfield(params.hrfParams, 'description')
         params.hrfParams.description = params.hrfModel;
      end
    
    
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
                    [scanParams, params] = getClassEVParamsGUI(thisView,params,defaultParams);
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
