function overlayList = selectOverlays(thisView,varargin)
% overlayList = selectOverlays(thisView,[title],[preselected]);
%
%   this function is deprecated, use overlayList = selectInList(thisView,'overlay',title,preselected)
%
% $Id$ 

overlayList = selectInList(thisView,'overlay',varargin);
mrWarnDlg('(selectOverlays) selectOverlays is deprecated. Please use ''selectInList(view,''overlays'',...)'' instead.');
