% makeEmptyROI.m
%
%        $Id:$ 
%      usage: makeEmptyROI(v,<scanNum>,<groupNum>)
%         by: justin gardner
%       date: 12/31/11
%    purpose: creates an empty roi with coordiantes set for teh scan and group
%
function roi = makeEmptyROI(v,varargin)

% check arguments
if nargin < 1
  help makeEmptyROI
  return
end

if ~isview(v),disp(sprintf('(makeEmptyROI) Invalid view passed in'));return,end

% get arguments
scanNum = [];groupNum = [];
getArgs(varargin,{'scanNum=[]','groupNum=[]'});

% make a name
if ieNotDefined('name')
  % go through roi names and get the largest numbered
  % roi name, i.e. ROI4 then make then new name ROI5
  maxnum = 0;
  roiNames = viewGet(v,'roiNames');
  for i = 1:length(roiNames)
    if regexp(roiNames{i},'^ROI\d+$')
      maxnum = max(maxnum,str2num(roiNames{i}(4:end)));
    end
  end
  name=sprintf('ROI%.0f',maxnum+1);
end
roi.name = name;
roi.voxelSize = viewGet(v,'scanVoxelSize',scanNum,groupNum);
if viewGet(v,'scanSformCode',scanNum,groupNum)
  roi.xform = viewGet(v,'scanSform',scanNum,groupNum);
else
  roi.xform = viewGet(v,'scanQform',scanNum,groupNum);
end

[tf roi] = isroi(roi);





