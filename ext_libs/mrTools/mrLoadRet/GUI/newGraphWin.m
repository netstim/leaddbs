function h = newGraphWin
%
% newGraphWin
%
%  $Id$
%
% Opens a new window and adds its handle to the global MLR.graphFigure.
% Sets closeRequestFcn to clean up properly (by calling
% closeGraphWin) when the window is closed.
%
% djh, 3/3/98
% djh, 9/2005 updated to MLR 4.0

mrGlobals
h=figure;
MLR.graphFigure = [h MLR.graphFigure];
set(h,'CloseRequestFcn','closeGraphWin');
figloc = mrGetFigLoc(['graphFigure' int2str(mlrGetFignum(h))]);
if ~isempty(figloc)
  %deal with multiple monitors
  [whichMonitor,figloc]=getMonitorNumber(figloc,getMonitorPositions);
  set(h,'Position',figloc);
end
