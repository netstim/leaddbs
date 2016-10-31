% mrParamsDialogSelectScans.m
%
%      usage: paramsInfo = mrParamsDialogSelectScans(v,groupNum,<paramsInfo>,<includeList>)
%         by: justin gardner
%       date: 08/17/07
%    purpose: add parameterInfo to select scans
%       e.g.:
%             v = newView;
%             paramsInfo = mrParamsDialogSelectScans(v,1);
%             params = mrParamsDialog(paramsInfo);
%             if ~isempty(params)
%               scanNums = find(params.include);
%             end
% 
%             if you want to only have one scan selected (say scan 3)
%             at start up:
%
%             paramsInfo = mrParamsDialogSelectScans(v,1,{},3);
%
function paramsInfo = mrParamsDialogSelectScans(v, groupNum,paramsInfo,includeList)

% check arguments
if ~any(nargin == [2 3 4])
  help mrParamsDialogSelectScans
  return
end

% get number of scans
nScans = viewGet(v,'nScans',groupNum);

% default arguments
if ieNotDefined('paramsInfo'),paramsInfo = {};end
initScan = 1;
if ieNotDefined('includeList')
  defaultSelected(1:nScans) = 1;
else
  defaultSelected(1:nScans) = 0;
  defaultSelected(includeList) = 1;
  if length(includeList) == 1
    initScan = min(includeList,nScans);
  end
end
% get info about scans
for i = 1:nScans
  description{i} = viewGet(v,'description',i,groupNum);
  filename{i} = viewGet(v,'tseriesFile',i,groupNum);
  includeScans{i} = defaultSelected(i);
end
paramsInfo{end+1} = {'scanNum',initScan,'round=1','incdec=[-1 1]',sprintf('minmax=[1 %i]',nScans),'Scan selector. Use this to choose which scans to include'};
paramsInfo{end+1} = {'include',includeScans,'type=checkbox','group=scanNum', 'Check this to include a particular scan, uncheck to skip.'};
includeNum = length(paramsInfo);
paramsInfo{end+1} = {'filename',filename,'group=scanNum','type=String','editable=0','Filename of scan'};
paramsInfo{end+1} = {'description',description,'group=scanNum','type=String','editable=0','Description of scan'};
paramsInfo{end+1} = {'selectAll',[],'type=pushbutton','buttonString=Select all scans','callback',@selectAll,'callbackArg',includeNum,'Click to include all scans'};

%%%%%%%%%%%%%%%%%%%
%%   selectAll   %%
%%%%%%%%%%%%%%%%%%%
function val = selectAll(includeNum)

global gParams
val = [];

% set the all values
gParams.varinfo{includeNum}.value = 1;
for i = 1:length(gParams.varinfo{includeNum}.allValues)
  gParams.varinfo{includeNum}.allValues{i} = 1;
end
% and update gui
set(gParams.ui.varentry{includeNum},'Value',1);
