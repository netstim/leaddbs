function view = newView(viewType)
%
%  v = newView
%
% djh, 6/2004

if nargin == 0
  viewType = 'Volume';
end

view = [];

% make sure paths are fixed to not to conflict with vista
mlrPath mrTools

% Define and initialize global variable MLR.
mrGlobals
if isempty(MLR.session)
  disp(sprintf('(newView) Could not open a view for directory %s',pwd));
  return
end

% check to make sure that we are not opening up MLR
% in a different home directory from where it started
if ~strcmp(pwd,MLR.homeDir) && isfile('mrSession.mat')
  if mlrIsRunning(0)
    mrWarnDlg(sprintf('(newView) MLR is open to a session in directory %s, but your current directory is %s. If you meant to open a new session for the current directory and not for the one in the MLR global, you will need to run mrQuit -this will close your mrLoadRet and deleteView all open views',MLR.homeDir,pwd));
  else
    % restart
    mrQuit;
    mrGlobals
    if isempty(MLR.session)
      disp(sprintf('(newView) Could not open a view for directory %s',pwd));
      return
    end
  end
end

viewNum = length(MLR.views) + 1;
view.viewNum = viewNum;
view.viewType = viewType;

% Initialize anat
view.baseVolumes = struct([]);
view.curBase = [];

% Initialize analysis list
view.analyses = {};
view.curAnalysis = [];

% Initialize ROIs
view.ROIs = struct([]);
view.curROI = [];
view.prevROIcoords = '';%this is a trick to detect that no coordinates have been put in this field, not even empty coordinates
view.showROIs = 'all';

% Initialize curGroup
view.curGroup = 1;
view.curScan = 1;

% Figure handle
view.figure = [];

% Coordinates corresponding to the slice currently displayed
view.curslice.baseCoords = [];
view.curslice.overlayCoords = [];

% Add the new view to the list of views
MLR.views{viewNum} = view;

% validate view (add any optional fields);
[tf view] = isview(view);

% add it back to the globals, in case it has changed
MLR.views{viewNum} = view;

% add the caches
MLR.caches{viewNum} = [];

% init the caches
view = viewSet(view,'roiCache','init');
view = viewSet(view,'overlayCache','init');
view = viewSet(view,'baseCache','init');


return;
