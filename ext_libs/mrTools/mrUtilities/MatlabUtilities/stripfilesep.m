% stripfilesep.m
%
%        $Id$
%      usage: pathname = stripfilesep(pathname)
%         by: justin gardner
%       date: 12/18/07
%    purpose: strip the file separator off a path (usually returned
%             by uigetfile)
%
function pathname = stripfilesep(pathname)

% check arguments
if ~any(nargin == [1])
  help stripfilesep
  return
end

if length(pathname) && (pathname(end)==filesep)
  pathname = pathname(1:end-1);
end
