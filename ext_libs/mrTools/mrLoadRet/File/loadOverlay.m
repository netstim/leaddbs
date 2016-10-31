function view = loadOverlay(view,filename,startPathStr,name)
%
%        $Id$
% view = loadOverlay(view,[filename],[startPathStr],[name])
%
% Loads an overlay (parameter map) array and adds it to view.overlays.
%
% If filename is not specified, prompts user to select a file. If filename
% is specified, it is assumed to be relative to startPathStr which is
% assumed (if not specified) to be the current group's data directory
% (e.g., view.subdir/Overlays/name.mat). Filename can be a string
% specifying an overlay file or it can be a cell array of filenames to load
% multiple overlays at once.
%
% name: optional replacement name(s) for the analysis when loaded. If
% pathStr is a cell array, then name must be a cell array of the same
% length.
%
% The file must contain a structure or structures, each with the following
% fields:
% - name: string
% - function: string function name by which it was computed
% - reconcileFunction: string function name that reconciles overlay data
%   and params with tseries files. See corAnalReconcileParams for an
%   example.
% - params: structure specify arguments to function
% - data: cell array of [y x z] arrays
% - date: specifies when it was computed, typically generated using
%         datestr(now) 
% - range: [min max] values
% - clip: [min max] to be displayed/thresholded
% - colormap: 256x3 array of RGB values
% - alpha: transparency value for alphaSlider
%
% In addition the 'name' field is set to the variable name to ensure
% consistency.
%
% djh, 1/9/98
% 6/2005, djh, update to mrLoadRet-4

mrGlobals

% Path to overlays
if ieNotDefined('startPathStr')
  startPathStr = viewGet(view,'overlayDir');
  if isempty(startPathStr) %if startPathStr is empty, that means there are no loaded analysis
     mrWarnDlg('(loadOverlay) Cannot load an overlay without an analysis');
     return
  end
end


% Complete pathStr
if ieNotDefined('filename')
    pathStr = mlrGetPathStrDialog(startPathStr,'Choose one or more overlays','*.mat','on');
else
	if iscell(filename)
		pathStr = cell(size(filename));
		for p=1:length(pathStr)
			pathStr{p} = fullfile(startPathStr,filename{p});
      pathStr{p} = [stripext(pathStr{p}),'.mat'];
		end
	else
		pathStr = {[stripext(fullfile(startPathStr,filename)),'.mat']};
	end
end

if ieNotDefined('name')
  name = [];
elseif ~iscell(name)
  name = {name};
end

if ~iscell(pathStr)
    pathStr = {pathStr};
end


% Loop through the filenames, load each file, and add each of them as a new
% overlay.
for p = 1:length(pathStr)
	if exist(pathStr{p},'file')
		h = mrMsgBox(['Loading overlay: ',pathStr{p},'. Please wait']);
		s = load(pathStr{p});
		varNames = fieldnames(s);
		overlay = s.(varNames{1});
    if isempty(name)
      overlay.name = varNames{1};
    else
      overlay.name = name{p};
    end
		% Add it to the view
		view = viewSet(view,'newOverlay',overlay);
		mrCloseDlg(h);
	else
		mrWarnDlg(['Overlay ',pathStr{p},' not found.']);
	end
end

return;

% Test/debug
view = loadOverlay(MLR.views{1},'co');
view = loadOverlay(MLR.views{1},{'co','amp','ph'});
view = loadOverlay(MLR.views{1});

