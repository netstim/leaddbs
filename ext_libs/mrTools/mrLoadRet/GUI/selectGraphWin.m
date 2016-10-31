function h = selectGraphWin(noClear,graphWindowPreference)
%
% selectGraphWin
%
%  $Id$
%
%  usage: h = selectGraphWin(noClear,graphWindowPreference)
%
% if  preference graphWindowPreference is set to replace, goes through the stack of figures and find first figure that is an MLR figure
% if no figure or preference graphWindowPreference is 'Make new', makes a new one.
% option noclear: replaces the current figure  but does not clear it (overrides graphWin preference)
%
% djh, 3/3/98
% djh, 9/2005 updated to MLR 4.0
% jlg 9/2006 added option not to clear window

mrGlobals;

% if there is no mrLoadRet running, just return a figure
if isempty(MLR.views)
  h = figure;
  return
end

if ieNotDefined('noClear')
  noClear = 0;
end
if ieNotDefined('graphWindowPreference')
  graphWindowPreference = mrGetPref('graphWindow');
end


if strcmpi(graphWindowPreference,'replace') || noClear
  figureHandles = get(0,'Children') ;
  % get current figure from MLR variable
  h = find(ismember(figureHandles,MLR.graphFigure),1,'first');
  if h
    h = figureHandles(h);
    %set(0,'CurrentFigure',h); this is redundant with figure(h)
    figure(h);
    % Clear the figure
    if (noClear == 0)
      clf
    end  
    return
  end
end

%if we're still here, that means there is no current figure that's part of MLR
%if (strcmp(mrGetPref('graphWindowPreference'),'Make new') && ~noClear) || (isempty(h) || h(1) == 0) || ~ishandle(h)  
h = newGraphWin;

return;
