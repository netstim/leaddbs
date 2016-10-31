% getLastDir(pathStr)
%
%      usage: getLastDir(pathStr,<levels>)
%         by: justin gardner
%       date: 04/05/07
%    purpose: Gets the last dir from a path. If levels is set then
%             it returns more directory levels (default for levels is 1).
%
function lastDir = getLastDir(pathStr,levels)

% check arguments
if ~any(nargin == [1 2])
  help getLastDir
  return
end

if ieNotDefined('levels'),levels = 1;end
if levels == 0,lastDir = '';return;end

% remove trailing fileseparator if it is there
if length(pathStr) && (pathStr(end) == filesep)
  pathStr = pathStr(1:end-1);
end

% get last dir
[pathStr lastDir ext] = fileparts(pathStr);

% paste back on extension
lastDir = [lastDir ext];

lastDir = fullfile(getLastDir(pathStr,levels-1),lastDir);
