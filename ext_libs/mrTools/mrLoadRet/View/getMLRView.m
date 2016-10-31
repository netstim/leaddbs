% getMLRView.m
%
%        $Id$
%      usage: getMLRView()
%         by: justin gardner
%       date: 12/21/07
%    purpose: returns the current view that is being used by
%             mrLoadRet (i.e. the GUI). This is a quick shortcut
%             to "steal" the view so that you can get information
%             from the command line about what is being display
%
%             the preferred usage of this is to substitute this
%             function call every place that requires a view:
%
%             viewGet(getMLRView,....
%             viewSet(getMLRView,...
%              
%             viewSet returns a modified view, but also updates
%             the global, so if you always reget the view from
%             the global using getMLRView, this should give
%             expected results.
%
function v = getMLRView()

% check arguments
v = [];
if ~any(nargin == [0])
  help getView
  return
end

mrGlobals;

% check for veiws
if isempty(MLR) || ~isfield(MLR,'views')
  disp(sprintf('(getView) No current MLR is running'));
  return;
end

% now go through each view and look for ones with valid figure numbers
viewsWithFigure = 0;
for i = 1:length(MLR.views)
  if isfield(MLR.views{i},'figure') && ~isempty(MLR.views{i}.figure)
    v = MLR.views{i};
    viewsWithFigure = viewsWithFigure + 1;
  end
end

if viewsWithFigure > 1
  disp(sprintf('(getView) Multiple figures open, returning last one'));
end

