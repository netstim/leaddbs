% mlrAnatDBGetUsername.m
%
%        $Id:$ 
%      usage: userName = mlrAnatDBGetUsername()
%         by: justin gardner
%       date: 07/05/15
%    purpose: Get username from mercurial
%
function userName = mlrAnatDBGetUsername()

% check arguments
if ~any(nargin == [0])
  help mlrAnatDBGetUsername
  return
end

% get user name
[status,userName] = system('hg config ui.username');
userName = strtrim(userName);


