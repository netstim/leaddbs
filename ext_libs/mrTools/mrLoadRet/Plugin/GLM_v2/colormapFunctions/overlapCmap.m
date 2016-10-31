function cmap = overlapCmap(numGrays,numColors)
%
% cmap = overlapCmap(numGrays,numColors)
% 
% Makes colormap array with:
%   gray scale - 1:numGrays
%   3-color overlap map - numGrays+1:numGrays+numColors
%
% ras 2/04

  if ~exist('numGrays','var')
    numGrays=128;
  end
  if ~exist('numColors','var')
    numColors=96;
  end

  cmap = zeros(numGrays+numColors,3);
  cmap(1:numGrays+numColors,:) = [gray(numGrays);zeros(numColors,3)];
  cmap(numGrays+1:numGrays+numColors/3,1) = .7;
  cmap(numGrays+numColors/3+1:numGrays+2*numColors/3,2) = .7;
  cmap(numGrays+2*numColors/3:numGrays+numColors,[3]) = .7;

  return
