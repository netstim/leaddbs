% getNumSlicesAtATime.m
%
%        $Id$
%      usage: getNumSlicesAtATime(numVolumes,dims)
%         by: justin gardner
%       date: 05/18/07
%    purpose: returns how many slices to process at a time,
%             dependent on the preference maxBlocksize
%             numVolumes is the number of volumes in the scan
%             dims is the dimensions of the scan
%            
%             The number of slices will always be an integer value
%             >=1. To get the unrounded value, check the second
%             return argument
%
%             A third return argument returns the number of
%             numRowsAtATime to do, in cases when memory constraints
%             require less than a single slice
%
%             A fourth return argument if accepted will recommend
%             precision to do computation at, i.e. double or
%             single. If you do not set this output argument then
%             the function will return recommended number of
%             slices/rows based on the precision passed in.
%    
function [numSlicesAtATime rawNumSlices numRowsAtATime precision] = getNumSlicesAtATime(numVolumes,dims,precision)

% check arguments
if ~any(nargin == [2 3])
  help getNumSlicesAtATime
  return
end

% get the precision (# of bytes)
if ieNotDefined('precision'),precision = 'double';end

switch (precision)
 case {'double'}
  bytesPerNum = 8;
 case{'single'}
  bytesPerNum = 4;
 otherwise
  bytesPerNum = 8;
end

% calculate slices at a time
maxBlocksize = mrGetPref('maxBlocksize');
if ieNotDefined('maxBlocksize')
  maxBlocksize = 250000000;
end
rawNumSlices = maxBlocksize/(bytesPerNum*numVolumes*prod(dims(1:2)));
numSlicesAtATime = max(1,floor(rawNumSlices));

% if we have less than one slice at a time, then we have to return
% numRowsAtATime as less than the full size
if (rawNumSlices < 1)
  % do at least two rows at a time (to avoid matrix getting
  % singleton dimension which causes squeeze to remove it)
  numRowsAtATime = max(round(dims(1)*rawNumSlices),2);
else
  numRowsAtATime = dims(1);
end

% check to see if we need to decrease precision, only do this
% if the calling function checks the precision output
if (nargout==4) && (numRowsAtATime == 2) && (strcmp(precision,'double'))
  % reset to single precision and recalculate
  [numSlicesAtATime rawNumSlices numRowsAtATime] = getNumSlicesAtATime(numVolumes,dims,'single');
  disp(sprintf('(getNumSlicesAtATime) Maximum blocksize of %i for dimensions: [%i %i %i %i] reached. Switching to single precision.',maxBlocksize,dims(1),dims(2),dims(3),numVolumes));
  precision = 'single';
end
