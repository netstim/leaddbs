% mrDefaultsFilename.m
%
%      usage: mrDefaultsFilename()
%         by: justin gardner
%       date: 09/21/07
%    purpose: gets the name of the defaults file. this
%             is used so that the user can override the
%             default location ~/.mrDefaults by setting
%             the preference
%
%             setpref('mrLoadRet','mrDefaultsFilename','~/.myMrDefaultsLoc');
%
function defaultsFilename = mrDefaultsFilename()

% check arguments
if ~any(nargin == [0])
  help mrDefaultsFilename
  return
end

if ispref('mrLoadRet','mrDefaultsFilename')
  defaultsFilename = getpref('mrLoadRet','mrDefaultsFilename');
else
  if ispc
    homeDir = [getenv('HOMEDRIVE') getenv('HOMEPATH')];
  else
    homeDir = getenv('HOME');
  end
  defaultsFilename = fullfile(homeDir,'.mrDefaults');
end
defaultsFilename = [defaultsFilename '.mat'];



