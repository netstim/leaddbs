% cellToCommaDelimited.m
%
%        $Id$
%      usage: str = cellToCommaDelimited(cellarray)
%         by: justin gardner
%       date: 10/05/07
%    purpose: converted a cell array of strings into a comma
%             delimited list
%
function out = cellToCommaDelimited(in)

% check arguments
if ~any(nargin == [1])
  help cellToCommaDelimited
  return
end
in = cellArray(in);
out = '';
if isempty(in)
  return
end
out = in{1};
for i = 2:length(in)
  out = sprintf('%s,%s',out,in{i});
end


