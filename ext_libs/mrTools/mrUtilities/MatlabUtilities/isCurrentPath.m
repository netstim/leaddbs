% isCurrentPath.m
%
%        $Id$
%      usage: tf = isCurrentPath(pathname)
%         by: justin gardner
%       date: 12/18/07
%    purpose: sees whether path is current path or not
%
function tf = isCurrentPath(pathname)

% check arguments
if ~any(nargin == [1])
  help isCurrentPath
  return
end

% make sure that there is no file separator at end of pathname
pathname = stripfilesep(pathname);

% cheng against current path
tf = strcmp(pwd,pathname);

