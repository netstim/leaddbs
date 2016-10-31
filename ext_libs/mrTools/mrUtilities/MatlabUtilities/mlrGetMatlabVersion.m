% mlrGetMatlabVersion.m
%
%        $Id:$ 
%      usage: versionNum = mlrGetMatlabVersion()
%         by: justin gardner
%       date: 04/02/15
%    purpose: get matlab version as a number
%
function versionNum = mlrGetMatlabVersion()

% check arguments
if ~any(nargin == [0])
  help mlrGetMatlabVersion
  return
end

% get version
v = version;

% parse the string (why, oh, why is it a string anyway?)
[majorVersion v] = strtok(v,'.');
[minorVersion v] = strtok(v,'.');

% convert to number
majorVersion = str2num(majorVersion);
minorVersion = str2num(minorVersion);

% and return as a single number
versionNum = majorVersion + minorVersion/10;

