% mrQuit.m
%
%        $Id$ 
%      usage: mrQuit(<saveMrLastView>,<v>)
%         by: justin gardner
%       date: 07/11/08
%    purpose: Function used to quit the mrLoadRet viewer. This is called from
%             the quit menu item, but if you call from the command line directly without 
%             any arguments it will close all open views and the viewer. Note that if you
%             have multiple viewers open on the data set, it will use the last view
%             opened to save settings in mrLastView:
%
%             mrQuit
%
%             Set saveMrLastView to 0 if you do not want to save a mrLastView. i.e.
%
%             mrQuit(0);
%
%             jb 2011/08/03 took some code out into mrSaveLastView.m so that it can be called by other functions
%
%
function retval = mrQuit(saveMrLastView,v)

% check arguments
if ~any(nargin == [0 1 2])
  help mrQuit
  return
end

mlrPath mrTools;
mrGlobals;

if ~ieNotDefined('v') && strcmp(class(v),'matlab.ui.eventdata.WindowCloseRequestData')
    disp('(mrQuit) v not passed in [CloseRequestData instead]')
    v = [];
end

% look for the open figure view
if ieNotDefined('v')
  v = [];
  % go through and look for the view with the figure
  if isfield(MLR,'views')
    for i = 1:length(MLR.views)
      if ~isempty(MLR.views{i}) && isfield(MLR.views{i},'figure') && ~isempty(MLR.views{i}.figure)
	if isempty(v)
	  v = MLR.views{i};
	else
	  fprintf('(mrQuit) Multiple open views found. Using last view opened for mrLastView')
	  v = MLR.views{i};
	end
      end
    end
  end
end

if (ieNotDefined('saveMrLastView') || mlrGetFignum(saveMrLastView) ) && ~isempty(v)
  mrSaveView(v);
end

if isfield(MLR,'views') && ~isempty(MLR.views)
  % close graph figure, remembering figure location
  if ~isempty(MLR.graphFigure)
    mrSetFigLoc('graphFigure',get(MLR.graphFigure,'Position'));
    close(MLR.graphFigure);
    MLR.graphFigure = [];
  end
  % close view figures
  % keep a local copy of everything since
  % the last view that is deleted will
  % clear the MLR global
  views = MLR.views;
  viewCount = 0;
  for viewNum = 1:length(views)
    view = views{viewNum};
    if isview(view)
      viewCount = viewCount+1;
      delete(view.figure); % deletes object, but apparently we still need
      closereq;            % a closereq to remove the window in r2014b
    end
  end
  if viewCount > 1
    fprintf('(mrQuit) Closing %i open views',viewCount);
  end
  drawnow
  disppercent(-inf,sprintf('(mrQuit) Saving %s',mrDefaultsFilename));
  saveMrDefaults;
  disppercent(inf);
else
  if ~isempty(v)
    closereq;
  end
end
clear global MLR

% revert paths to what they were from before running mrTools
mlrPath revert


