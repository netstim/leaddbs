% mlrMakeAnalysis.m
%
%        $Id:$ 
%      usage: a = mlrMakeAnalysis(v,analysisName)
%         by: justin gardner
%       date: 11/10/14
%    purpose: Makes a basic analysis for use with MLR. This is a very
%             simple way to put your overlay into MLR (see also mrDispOverlay)
%             For example, say you have computed an overaly with values from
%             0 to 1 for each voxel in the current scan (i.e. dimensions
%             [x y s]. Then you can do:
%
%             % Get the current viewers v
%             v = getMLRView;
%
%             % create a random overlay for this view
%             overlay = rand(viewGet(v,'scanDims'));
%
%             % Make the most basic analysis to show this overlay
%             a = mlrMakeAnalysis(v,'myAnalysis',overlay);
%
%             % set it
%             v = viewSet(v,'newAnalysis',a);
%
%             % and view it
%             refreshMLRDisplay(v);
% 
%             and your overlay should appear in the viewer.
%
%             Note that this can be called iteratively to add
%             multiple overlays:
%
%             a = mlrMakeAnalysis(v,'myAnalysis',overlay);
%             a = mlrMakeAnalysis(v,a,overlay2,'overlay2name');
%    
%             You can also set various parameters of the overlays like
%             the colormap and the name of the overlay:
%             a = mlrMakeAnalysis(v,'myAnalysisName',overlay,'overlayName','myOverlayName','overlayColormap',hot(256));
%
%
function a = mlrMakeAnalysis(v,analysisName,overlay,varargin)

% check arguments
if nargin < 2
  help mlrMakeAnalysis
  return
end

% if analysisName is a string then create the analysis
if isstr(analysisName)
  a.groupName = viewGet(v,'groupName');
  a.name = analysisName;
  [tf a] = isanalysis(a);
else
  [tf a] = isanalysis(analysisName);
end

% check for valid analysis
if ~tf
  disp(sprintf('(mlrMakeAnalysis) Analysis is not valid'));
  return
end

% get the arguments
getArgs(varargin,{'overlayParams',[],'overlayRange',[],'overlayFunction',[],'overlayColormap',[],'overlayName',a.name,'overlayClip',[]});

% if there is an overlay then add that
if nargin >= 3
  % get scanNum
  scanNum = viewGet(v,'curScan');
  numScans = viewGet(v,'numScans');
  % validate overlay
  scanDims = viewGet(v,'scanDims');
  if ~isequal(scanDims,size(overlay))
    disp(sprintf('(mlrMakeAnalysis) Overlay dimensions [%s] should match scan dimensiosn [%s]',mlrnum2str(size(overlay),'sigfigs=0'),mlrnum2str(scanDims,'sigfigs=0')));
    return
  end
  % set the data
  for iScan = 1:numScans
    if iScan == scanNum
      o.data{iScan} = overlay;
    else
      o.data{iScan} = [];
    end
  end
  % set overlayFunction
  o.function = overlayFunction;
  % set group to be current group
  o.groupName = viewGet(v,'groupName');
  % set overlay name
  o.name = overlayName;
  % set this overlay to be for this scan
  if isempty(overlayParams) || ~isfield('scanList',overlayParams)
    overlayParams.scanNum = viewGet(v,'curScan');
  end
  % set params
  o.params = overlayParams;
  % set colormap
  if ~isempty(overlayColormap)
    o.colormap = overlayColormap;
  end
  % set overlayRange
  if isempty(overlayRange)
    overlayRange = [min(overlay(:)) max(overlay(:))];
  end
  o.range = overlayRange;
  if ~isempty(overlayClip)
    o.clip = overlayClip;
  end
  % validate overlay
  [tf o] = isoverlay(o);
  if ~tf
    disp(sprintf('(mlrMakeAnalysis) Could not validate overlay'));
  else
    % add to analysis
    if ~isfield(a,'overlays') || isempty(a.overlays)
      a.overlays = o;
    else
      a.overlays(end+1) = o;
    end
  end
  a.curOverlay = 1;
end
  

