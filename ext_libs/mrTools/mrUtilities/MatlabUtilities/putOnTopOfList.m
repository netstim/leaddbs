% putOnTopOfList.m
%
%      usage: putOnTopOfList(topVal,list)
%         by: justin gardner
%       date: 04/03/07
%    purpose: puts the topVal on top of the list
%
% e.g.
%putOnTopOfList('topItem',{'oneItem','twoItem','topItem','threeItem'});
%
function outList = putOnTopOfList(topVal,inList)

% check arguments
if ~any(nargin == [2])
  help putOnTopOfList
  return
end

% handle empty lists correctly
if isempty(inList)
  inList = {};
end
if isempty(topVal)
  outList = inList;
  return
end

% make sure inList is a cell array
inList = cellArray(inList);

% put the top val at the top of list
outList{1} = topVal;

% and add everybody else
for i = 1:length(inList)
  if ~isequalwithequalnans(inList{i},topVal)
    outList{end+1} = inList{i};
  end
end


