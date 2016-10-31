% cellArray.m
%
%      usage: var = cellArray(var,<numLevels>)
%         by: justin gardner
%       date: 04/05/07
%    purpose: when passed a single structure it returns it as a
%    cell array of length 1, if var is already a cell array just
%    passes it back. numLevels controls how many levels of
%    cell array you want, usually this would be one, but
%    if you wanted to make sure you have a cell array of cell
%    arrays then set it to two.
%  
% 
% e.g.:
% c = cellArray('this')
%
function var = cellArray(var,numLevels)

% check arguments
if ~any(nargin == [1 2])
  help cellArray
  return
end

if ieNotDefined('numLevels'),numLevels = 1;,end

% deal with case in which we are trying to make a multi
% level cell array into a single cell array
if numLevels == -1
  if iscell(var)
    tmp = {};
    for i = 1:length(var)
      tmpi = cellArray(var{i},-1);
      tmp = {tmp{:} tmpi{:}};
    end
    var = tmp;
  else
    tmp{1} = var;
    var = tmp;
  end
  return
end
  
% make the variable name
varName = 'var';

for i = 1:numLevels
  % for each level make sure we have a cell array
  % if not, make it into a cell array
  if ~iscell(eval(varName))
    tmp = var;
    clear var;
    var{1} = tmp;
  end
  % test the next level
  varName = sprintf('%s{1}',varName);
end



