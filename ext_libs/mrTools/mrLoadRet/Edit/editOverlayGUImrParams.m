% editOverlayGUImrParams.m
%
%        $Id$	
%      usage: editOverlayGUImrParams(viewNum)
%         by: eli merriam
%       date: 09/04/07
%    purpose: 
%
function retval = editOverlayGUImrParams(viewNum)

% check arguments
  if ~any(nargin == [1])
    help editOverlayGuimrParams
    return
  end

  v = viewGet(viewNum,'view');

  % Get the original overlay
  analysisNum = viewGet(v,'currentAnalysis');
  if isempty(analysisNum),mrWarnDlg('(editOverlayGUI) No current analysis');return,end
  overlayNum = viewGet(v,'currentOverlay', analysisNum);
  if isempty(overlayNum),mrWarnDlg('(editOverlayGUI) No current overlay');return,end
  overlayRange = viewGet(v,'overlayRange', overlayNum, analysisNum);
  overlayClip = viewGet(v,'overlayClip', overlayNum, analysisNum);
  overlayName = viewGet(v, 'overlayName', overlayNum, analysisNum);
  alphaOverlay = viewGet(v,'alphaOverlay');
  alphaOverlayExponent = viewGet(v,'alphaOverlayExponent');
  interrogator = viewGet(v,'interrogator',overlayNum,analysisNum);
  
  % colormaps
  colormaps = {'default','hot','hsv','pink','cool','bone','copper','flag','gray','grayCirc','twoCondCmap','twoCondCircCmap','hsvDoubleCmap','cmapExtendedHSV','overlapCmap','redGreenCmap','rygbCmap','bicolorCmap' 'coolCmap','hotColdCmap'};
  altColormaps = viewGet(v,'colormaps');
  if ~isempty(altColormaps)
    colormaps = {colormaps{:} altColormaps{:}};
  end
  
  % set up params dialog
  paramsInfo = {};
  paramsInfo{end+1} = {'overlayCmap', colormaps,'type=popupmenu','List of possible colormaps','callback',@mrCmapCallback,'callbackArg',v};
  paramsInfo{end+1} = {'userDefinedCmap','','Allows you to call a user defined function to set the overla colormap','callback',@mrCmapCallback,'callbackArg',v};
  paramsInfo{end+1} = {'numColors', 256, 'first argument to the colormap function','callback',@mrCmapCallback,'callbackArg',v};
  paramsInfo{end+1} = {'numGrays', 0, 'second argument to the colormap function','callback',@mrCmapCallback,'callbackArg',v};
  paramsInfo{end+1} = {'flipCmap', 0, 'type=checkbox', 'check this box to reverse the direction of the colormap','callback',@mrCmapCallback,'callbackArg',v};
  paramsInfo{end+1} = {'shiftCmap', 0, 'incdec=[-16 16]', 'shift the colormap -- this can be useful for retinotopy scans with circular colormaps','callback',@mrCmapCallback,'callbackArg',v}; 
  paramsInfo{end+1} = {'overlayCtype', {'normal', 'setRangeToMax', 'setRangeToMaxAroundZero'}, 'type=popupmenu','setRangeToMax scales the colormap to overlayMin-overlayMax, as in R2 maps','callback',@mrCmapCallback,'callbackArg',v};
  paramsInfo{end+1} = {'overlayRange', overlayRange, 'The lower and upper bound on the colormap','callback',@mrCmapCallback,'callbackArg',v};
  paramsInfo{end+1} = {'overlayClip', overlayClip, 'The lower and upper clip points on the colormap','callback',@mrCmapCallback,'callbackArg',v};
  paramsInfo{end+1} = {'interrogator', interrogator, 'Set the interrogator function name','callback',@mrCmapCallback,'callbackArg',v};
  paramsInfo{end+1} = {'alphaOverlay', alphaOverlay, 'You can specify the name of another overlay in the analysis to use as an alpha map. For instance, you might want to display one overlay with the alpha set to the r2 or the coherence value.','callback',@mrCmapCallback,'callbackArg',v};
 paramsInfo{end+1} = {'alphaOverlayExponent', alphaOverlayExponent,'incdec=[-0.1 0.1]','If you are using an alphaOverlay, this sets an exponent on the alphaOverlay to pass through. For example, if you just want the value from the overlay to be the alpha value then set this to 1. If you want to have it so that lower values get accentuated (usually this is the case), set the exponent to less than 1, but greater than 0. The alpha values are passed through the function alpha = alpha.^alphaOverlayExponent','callback',@mrCmapCallback,'callbackArg',v};
  paramsInfo{end+1} = {'setAll', 0, 'type=pushbutton','Set all overlays to have them same settings','callback',@mrCmapSetAllCallback,'callbackArg',v,'buttonString=Set all overlays','passParams=1'};
%   paramsInfo{end+1} = {'overlayName', overlayName, 'The name for the overlay (e.g., co, am, or ph)'};

  % display dialog
  mrParamsDialog(paramsInfo,'Change overlay colormap');

  return;


function mrCmapCallback(v,params)

  disppercent(-inf,'(editOverlayGUImrParams) Recomputing overlay');
  v = viewSet(v,'overlayCache','init');
  
  % get the current overlay
  analysisNum = viewGet(v,'currentAnalysis');
  overlayNum = viewGet(v,'currentOverlay',analysisNum);
  o = viewGet(v, 'overlay', overlayNum, analysisNum);

  % set which color cmap to use
  if ~strcmp(params.overlayCmap, 'default')
    if sum(strcmp(params.overlayCmap, {'hsvDoubleCmap','cmapExtendedHSV','cmapHSV','overlapCmap','redGreenCmap','rygbCmap','bicolorCmap','coolCmap'}))
      o.colormap = eval(sprintf('%s(%i,%i)', params.overlayCmap, params.numGrays, params.numColors));
    else
      o.colormap = eval(sprintf('%s(%i)', params.overlayCmap, params.numColors));
    end
  end

  % see if we need to call a function
  if ~isempty(params.userDefinedCmap)
    % look for the m function
    if exist(sprintf('%s.m',params.userDefinedCmap))
      colormap = eval(sprintf('%s(%i)',params.userDefinedCmap,params.numColors));
      if isequal(size(colormap),[params.numColors 3])
	o.colormap = colormap;
      else
	disp(sprintf('(editOverlay) Function %s must return a %ix%i array',params.userDefinedCmap,params.numColors,3));
      end
    end
  end
    
  % flip the cmap
  if params.flipCmap
    o.colormap = flipud(o.colormap);
  end

  % shift the cmap
  if params.shiftCmap
    o.colormap = circshift(o.colormap, params.shiftCmap);
  end

  % scale to max, or not
  o.colormapType = params.overlayCtype;

  % set the overlay range & clip
  o.range = [params.overlayRange(1) params.overlayRange(2)];
  o.clip = [params.overlayClip(1) params.overlayClip(2)];
  

%   % set the name of the overlay
%   o.name = params.overlayName;
  
  o.alphaOverlay = params.alphaOverlay;
  o.alphaOverlayExponent = params.alphaOverlayExponent;

  o.interrogator = params.interrogator;
  
  % set the new overlay
  v = viewSet(v,'newOverlay', o);

  % and refresh
  refreshMLRDisplay(v.viewNum);

  disppercent(inf);
  return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   mrCmapSetAllCallback   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = mrCmapSetAllCallback(v,params)
retval = [];
hWaitBar = mrWaitBar(0,'(editOverlayGUImrParams) Setting overlays to have same parameters');

% get the current overlay
currentOverlayNum = viewGet(v,'curOverlay');
currentOverlayName = viewGet(v,'overlayname',currentOverlayNum);
currentOverlay = viewGet(v,'overlay',currentOverlayNum);
nOverlays = viewGet(v,'numOverlays');

% a list of fields which should be copied
copyFields = {'alpha','alphaOverlayExponent','interrogator','clip','colormap','colormapType','range'};
% now go through each overlay and set the fields the same one

for onum = 1:nOverlays
  o = viewGet(v,'overlay',onum);
  % copy fields from the current overlay
  for fnum = 1:length(copyFields)
    o.(copyFields{fnum}) = currentOverlay.(copyFields{fnum});
  end
  % now set it back
  v = viewSet(v,'newOverlay',o);
  mrWaitBar(onum/nOverlays,hWaitBar);
end

% set back the current overlay to the original one
v = viewSet(v,'curOverlay',viewGet(v,'overlayNum',currentOverlayName));

mrCloseDlg(hWaitBar);

%%%%%%%%%%%%%%%%%%%%%%%%%
%%   mrEpiMovieClose   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function mrCmapClose

  return;
  

  
  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% list of useful cmaps taken from MLR3 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cmap = grayCirc(numGrays)
  if ~iseven(numGrays);
    numGrays = numGrays + 1;
  end
  
  cmapA = gray(numGrays/2 + 1);
  cmapB = flipud(cmapA);
  
  cmap = cat(1,cmapA(1:numGrays/2,:), cmapB(1:numGrays/2,:));

function cmap = hsvDoubleCmap(numGrays,numColors,symmetric)
%
% cmap = hsvDoubleCmap(numGrays,numColors,symmetric)
% 
% Makes colormap array with:
%   gray scale - 1:numGrays
%   hsv colors - numGrays+1:numGrays+numColors
%
% Double wrapping the colors is useful for visualize retinotopy
% so that 180 deg from one hemifield maps onto a full hsv
% colormap (instead of just half of it).  If symmetric, flip
% second hsv cmap
%
% djh 1/98

  if ~exist('numGrays','var')
    numGrays=128;
  end
  if ~exist('numColors','var')
    numColors=96;
  end
  if ~exist('symmetric','var')
    symmetric=1;
  end

  cmap = zeros(numGrays+numColors,3);
  if symmetric
    cmap(1:numGrays+numColors,:) = ...
        [gray(numGrays);
         hsv(floor(numColors/2));
         flipud(hsv(ceil(numColors/2)))];
  else
    cmap(1:numGrays+numColors,:) = ...
        [gray(numGrays);
         hsv(floor(numColors/2));
         hsv(ceil(numColors/2))];
  end


  
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

function cmap = redGreenCmap(numGrays,numColors)
%
% cmap = redGreenCmap(numGrays,numColors)
% 
% Makes colormap array with:
%   gray scale - 1:numGrays
%   redGreen ramp - numGrays+1:numGrays+numColors
%
% djh 1/98

  if ~exist('numGrays','var')
    numGrays=128;
  end
  if ~exist('numColors','var')
    numColors=96;
  end

  cmap = zeros(numGrays+numColors,3);
  cmap(1:numGrays,:) = gray(numGrays);
  for i=numGrays+1:numGrays+numColors
    cmap(i,:) = ...
        [((i-numGrays)/numColors)^.5, (1-(i-numGrays)/numColors)^.5, 0];
  end


function cmap = rygbCmap(numGrays,numColors)
%
% cmap = rygbCmap(numGrays)
% 
% Makes colormap array with:
%   gray scale - 1:numGrays
%   red, yellow, green, blue - next 4 colors
%
% djh 1/98

  if ~exist('numGrays','var')
    numGrays=128;
  end

  nRs=floor(1/4*numColors);
  nYs=floor(1/2*numColors)-nRs;
  nGs=floor(3/4*numColors)-nRs-nYs;
  nBs=numColors-nGs-nYs-nRs;

  cmap = zeros(numGrays+numColors,3);
  cmap(1:numGrays,:) = gray(numGrays);
  cmap(numGrays+1:numGrays+nRs,:) = ones(nRs,1)*[1 0 0];
  cmap(numGrays+nRs+1:numGrays+nRs+nYs,:) = ones(nYs,1)*[1 1 0];
  cmap(numGrays+nRs+nYs+1:numGrays+nRs+nYs+nGs,:) = ones(nGs,1)*[0 1 0];
  cmap(numGrays+nRs+nYs+nGs+1:numGrays+numColors,:) = ones(nBs,1)*[0 0 1];


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

  
function cmap = twoCondCmap(numGrays)
  minx = 0;maxx = 0.5;
  x = minx:(maxx-minx)/(numGrays/2):maxx;
  y = mygauss([1 0 0.25 0],x);
  
  red = [zeros(length(y),1) y' y'];
  y = fliplr(y);
  blue = [y' y' zeros(length(y),1)];
  
  cmap = 1-[red;blue];
return

function cmap = hotColdCmap(n)


  h = hot(floor(n/2));
  c(:,1) = h(:,3);
  c(:,2) = h(:,2);
  c(:,3) = h(:,1);
  
  if iseven(h)
    cmap = [flipud(c);h];
  else
    cmap = [flipud(c);[0 0 0];h];
  end

return

function cmap = twoCondCircCmap(numGrays)
  if isodd(numGrays)
    numGrays = numGrays +1;
  end
  cmapA = twoCondCmap(numGrays/2);
  cmapB = flipud(cmapA);
  cmap = cat(1,cmapA, cmapB);
return
  
function cmap = cmapExtendedHSV(numGrays,numColors,range)
%
% cmap = cmapExtendedHSV([numGrays=128],[numColors=96],[range=query user])
% 
% Makes a special purpose colormap that is useful for representing rotating
% wedge data.  The map created here are hsv maps where all of the colors
% can be placed in a subsection of the full color map.  In this way, the
% full range of colors spans less than 2pi, like in the double color map.
% Rather than compressing by a complete factor of 2, like the double color
% map, the compression factor can be a bit smaller.
%
%   There are numGrays gray scale entries.  They occupy the first part of the cmap, 1:numGrays
%   There are numColors hsv colors.  They fill the map entries following the gray, 
%   from (numGrays+1):numGrays+numColors
%   When the range is 1, the hsv map is hsv(numColors) and we insert it
%   into the cmap.
%   When the range is 1.5, we compute tmp = hsv(numColors/1.5) and we
%   create [tmp,tmp(1:needed)] to fill up numColors entries.
%
% Examples:
%   cmap = cmapExtendedHSV(128,96,1.1);
%   cmap = cmapExtendedHSV(128,96);    -- 128 gray levels, 96 color levels,
%          query use for compression
%   cmap = cmapExtendedHSV;
%   cmap = cmapExtendedHSV(128,96,2);  -- Same as hsvDoubleCmap
%


  if ~exist('numGrays','var'),    numGrays=128; end
  if ~exist('numColors','var'),   numColors=96; end
  if ~exist('range','var'),       range = readRange;  end

  if (range < 1) | (range > 2),  error('Range must be betweem 1 and 2.'); end

  hsvColors = round(numColors/range);
  hsvColorsExtra = numColors - hsvColors;
  hsvMap = cmapHSV(hsvColors);

  cmap = zeros(numGrays+numColors,3);

  % If you want the map symmetric at the boundary, you should do this.
  % We could trap range == 2 and do it then ... which would be backwards
  % compatible?
  % cmap = [gray(numGrays); hsvMap; flipud(hsvMap(1:hsvColorsExtra,:))];
  cmap = [hsvMap; hsvMap(1:hsvColorsExtra,:)];
  shiftSize = round(hsvColorsExtra/2);
  hsvMap = circshift(cmap,shiftSize);

  cmap = [gray(numGrays); cmap];

  return;
  
function hsvMap = cmapHSV(hsvColors);
%
%   hsvMap = cmapHSV(hsvColors);
% 
%Author: AB/BW
%Purpose:
%   Create an hsv map such that the magenta/blue boundary is in the middle.
%   This makes manipulation of the colors easier for wedge maps.
%

% shiftSize to make magenta  the middle color is in the center.
  magenta = [.1 0 1];
  mp = hsv(hsvColors);

  % Obscure bit of code.  Have fun.
  [val,idx] = min(sum( abs(mp - repmat(magenta,hsvColors,1))'));
  shiftSize = -round((idx - hsvColors/2));

  hsvMap = circshift(mp,shiftSize);

  return;

function cmap = coolCmap(numGrays,numColors)
%
% cmap = coolCmap(numGrays,numColors)
% 
% Makes colormap array with:
%   gray scale - 1:numGrays
%   cool colors - numGrays+1:numGrays+numColors
%
% djh 1/98

  if ~exist('numGrays','var')
    numGrays=128;
  end
  if ~exist('numColors','var')
    numColors=96;
  end

  cmap = zeros(numGrays+numColors,3);
  cmap(1:numGrays+numColors,:) = [gray(numGrays);cool(numColors)];



%----------------------------------------
function range = readRange

  prompt={'Enter compression range for the hsv map (2 = double color map)'};
  def={'1.2'};
  dlgTitle='Color map compression factor';
  lineNo=1;
  range=inputdlg(prompt,dlgTitle,lineNo,def);
  range = str2num(range{1});

  return;


% Gaussian
%
% usage: gauss(p,X,Y);
%   p is an array of parameters:
%     p(1) = height of Gaussian
%     p(2) = center x
%     p(3) = center y
%     p(4) = SD in x dimension
%     p(5) = SD in y dimension
%     p(6) = offset
%     p(7) = rotation in radians
%   X and Y are the position on which to evaluate the gaussian
%     to evaluate at a matrix of points use,
%     e.g. [X,Y] = meshgrid(-1:.1:1,-1:.1:1);
%
%  the function can also be called as follows for 1D
%  usage: gauss(p,X);
%
%     p(1) = height of Gaussian
%     p(2) = center
%     p(3) = SD
%     p(4) = offset
% 
%   by: justin gardner
% date: 6/6/97
function G=mygauss(p,X,Y)

% 2D Gaussian 
if nargin == 3

  % rotate coordinates
  % note that the negative sign is because
  % we are rotating the coordinates and not the function
  X1 = cos(-p(7)).*X - sin(-p(7)).*Y;
  Y1 = sin(-p(7)).*X + cos(-p(7)).*Y;

  % calculate the Gaussian
  G = p(1) * exp(-((((X1-p(2)).^2)/(2*p(4)^2))+(((Y1-p(3)).^2)/(2*p(5)^2))))+p(6);
  
% 1D Gaussian
elseif nargin == 2

  % calculate the Gaussian
  G = p(1) * exp(-(((X-p(2)).^2)/(2*p(3)^2)))+p(4);
  
else 
   % usage error
   disp('USAGE: gauss(parameters, X, Y)');
end
