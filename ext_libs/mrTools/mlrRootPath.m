% mlrRootPath.m
%
%        $Id:$ 
%      usage: mlrRootPath()
%         by: justin gardner
%       date: 11/20/14
%    purpose: returns the root path for mlr
%
function retval = mlrRootPath()

% check arguments
if ~any(nargin == [0])
  help mlrRootPath
  return
end

retval = fileparts(which('mlrRootPath'));

