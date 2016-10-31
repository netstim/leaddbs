% getext.m
%
%      usage: getext(filename)
%         by: justin gardner
%       date: 02/08/05
%    purpose: returns extension if it exists
%
function retval = getext(filename,separator)

if ~any(nargin == [1 2])
  help getext;
  return
end

if ~exist('separator'),separator = '.';,end
retval = '';
dotloc = findstr(filename,separator);
if (length(dotloc) > 0) && (dotloc(length(dotloc)) ~= length(filename))
  retval = filename(dotloc(length(dotloc))+1:length(filename));
end
