function cmap = bicolorCmap(numGrays,numColors)
%
% cmap = bicolorCmap(numGrays,numColors)
% 
% Makes colormap array with:
%   gray scale - 1:numGrays
%   cool colors - values in which map < 0
%   black - values in which map=0
%   hot colors - values in which map > 0
%
% This is useful for plotting contrast maps in which both
% positive and negative effects are displayed (for related updates, see
% loadParameterMap, computeContrastMap).
%
% djh 1/98
% ras, 03/04, off of hotCmap
%   numGrays  = 64;
%   numColors = 64;

%   hi = max(max(max(map)));
%   lo = min(min(min(map)));

  
% ATTN:  for now just assume that range is [-1 1]
  hi = 1;
  lo = -1;
  
  rng = linspace(lo,hi,numColors);
  
  %   if lo > 0 % all positive
  %     colors = hot(numColors);
  %   elseif hi < 0 % all negative 
  %     colors = cool(numColors);
  %   else        % crosses zero
  colors = zeros(numColors,3);
  neg = length(find(rng < 0));
  colors(neg,:) = [0 0 0]; % black when crosses
  colors(1:neg-1,:) = flipud(cool(neg-1));
  colors(neg+1:end,:) = hot(numColors-neg);
  %   end
  
  cmap = zeros(numGrays+numColors,3);
  cmap(1:numGrays+numColors,:) = [gray(numGrays);colors];
  
  return

  
