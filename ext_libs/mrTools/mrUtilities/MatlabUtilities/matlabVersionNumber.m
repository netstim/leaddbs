% matlabVersionNumber.m
%
%        $Id:$ 
%      usage: [majorVersion minorVersion] = matlabVersionNumber()
%         by: justin gardner
%       date: 05/08/13
%    purpose: Get matlab version number
%
function [majorVersion minorVersion] = matlabVersionNumber()

% check arguments
if ~any(nargin == [0])
  help matlabVersionNumber
  return
end

majorVersion = [];
minorVersion = [];

% get matlab version string
matlabVersion = ver('Matlab');

% read the string (note that it can be like 7.1.1 which doesn't process well through str2num)
versionNumbers = textscan(matlabVersion.Version,'%f.%f');

% if we have the first number it is the major version
if length(versionNumbers) >= 1
  majorVersion = versionNumbers{1};
end

% if we have the second number it is the minor version
if length(versionNumbers) >= 2
  minorVersion = versionNumbers{2};
end

% if minor version is empty set it to 0
if isempty(minorVersion)
  minorVersion = 0;
end
  


