% projectOutMeanVectorParams.m
%
%        $Id$
%      usage: projectOutMeanVectorParams(v,<defaultParams>)
%         by: justin gardner
%       date: 05/02/08
%    purpose: Get the parameters for projectOutMeanVectorParams
%             Set defaultParams=1 if you just want to get default params
%
function params = projectOutMeanVectorParams(v,defaultParams)

% check arguments
if ~any(nargin == [1 2])
  help projectOutMeanVectorParams
  return
end

params = [];

if ~ieNotDefined('defaultParams') && (defaultParams == 1)
  params.sourceROI = [];
  params.targetROI = [];
  return
end

% ask the user which rois to use as the source ROI
roiNames = viewGet(v,'roiNames');
roiNames = putOnTopOfList('Whole scan',roiNames);
source = buttondlg('Select source roi(s) from which mean vector will be computed',roiNames);
if sum(source)==0,return,end

% ask about target ROI
target = buttondlg('Select target roi(s) from which mean vector will be projected out',roiNames);
if sum(target)==0,return,end

if source(1)
  % if user selected whole scan then sourceROI should be set as empty
  params.sourceROI = [];
else
  % otherwise get a cell array of the roi names for the source
  params.sourceROI = {};
  for i = 1:length(source)
    if source(i)
      params.sourceROI{end+1} = roiNames{i};
    end
  end
end

if target(1)
  % if user selected whole scan then sourceROI should be set as empty
  params.targetROI = [];
else
  % otherwise get a cell array of the roi names for the source
  params.targetROI = {};
  for i = 1:length(target)
    if target(i)
      params.targetROI{end+1} = roiNames{i};
    end
  end
end


