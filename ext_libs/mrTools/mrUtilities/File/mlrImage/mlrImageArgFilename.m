% mlrImageArgFilename.m
%
%      usage: filename = mlrImageArgFilename(imageArg)
%         by: justin gardner
%       date: 09/10/11
%    purpose: Takes an imageArg returned by mlrImageParseArgs and
%             retruns its filename - or some other string that
%             identifies it - this is useful in programs that use
%             mlrImageParseArgs to report an error opening file
%
function filename = mlrImageArgFilename(imageArg)

filename = [];

% check arguments
if ~any(nargin == [1])
  help mlrImageArgFilename
  return
end

if isstr(imageArg)
  filename = imageArg;
elseif isstruct(imageArg) 
  if isfield(imageArg,'filename')
    filename = imageArg.filename;
  elseif isfield(imageArg,'data')
    filename = 'passed in data';
  end
end



