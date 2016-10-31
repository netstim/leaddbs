% reshapeTSeries
%
%      usage: tseries = reshapeTSeries(tseries)
%         by: eli merriam
%       date: 03/22/07
%    purpose: input is an (X,Y,T) slice, output is a reshaped (X*Y,T) matrix
%	$Id$	
%

function tseries = reshapeTSeries(tseries)

% check arguments
if ~any(nargin == [1])
  help reshapeTSeries
  return
end

% Reformat to nFrames x nVoxels
tseries = squeeze(tseries);
dims = size(tseries);
nFrames = dims(3);
sliceDims = dims(1:2);
nVoxels = prod(sliceDims);
tseries = reshape(tseries,[nVoxels,nFrames])';
