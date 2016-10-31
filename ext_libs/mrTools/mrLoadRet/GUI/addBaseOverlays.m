% addBaseOverlays.m
%
%        $Id:$ 
%      usage: overlays = addBaseOverlays(v,baseNum,overlays)
%         by: justin gardner
%       date: 11/09/14
%    purpose: add base overlays (those are overlays that are computed
%             in base voxels - like for fascicle stuff) to overlays -
%             should be called after computeoverlays (see refreshMLRDisplay)
%
function overlays = addBaseOverlays(v,baseNum,overlays)

% check arguments
if ~any(nargin == [3])
  help addBaseOverlays
  return
end

% get baseOverlay
baseOverlay = viewGet(v,'baseOverlay',baseNum);
if isempty(baseOverlay),return,end

% get size of base
baseDims = viewGet(v,'baseDims',baseNum);

% get baseOverlayAlpha
baseOverlayAlpha = viewGet(v,'baseOverlayAlpha',baseNum);
if isscalar(baseOverlayAlpha)
  baseOverlayAlpha = ones([baseDims(1:2) 3])*baseOverlayAlpha;
end

% figure out what overlay number to add to
if isempty(overlays.RGB)
  overlayNum = 1;
else
  overlayNum = size(overlays.RGB,4)+1;
end

% maken the baseOverlay in voxelsxRBG format
if isstr(baseOverlay)
  baseOverlay = color2RGB(baseOverlay);
end

% if it is a [R G B] then just repeat for every voxedl
if length(baseOverlay) == 3
  overlays.RGB(1:baseDims(1),1:baseDims(2),1,overlayNum) = baseOverlay(1);
  overlays.RGB(1:baseDims(1),1:baseDims(2),2,overlayNum) = baseOverlay(2);
  overlays.RGB(1:baseDims(1),1:baseDims(2),3,overlayNum) = baseOverlay(3);
else
  if ~isequal(size(baseOverlay),[baseDims(1:2) 3])
    disp(sprintf('(addBaseOverlays) Expecting baseOverlay to be XxYx3'));
  else
    overlays.RGB(:,:,:,overlayNum) = baseOverlay;
  end
end

% add alpha map
overlays.alphaMaps(:,:,:,overlayNum) = baseOverlayAlpha;

% add bogus cmap (so that refreshMLRDisplay has something)
overlays.cmap = [];
overlays.range = [0 1];

