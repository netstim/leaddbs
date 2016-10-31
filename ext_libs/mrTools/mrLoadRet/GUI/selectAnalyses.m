function analysesList = selectAnalyses(thisView,varargin)
% analysesList = selectAnalyses(thisView,[title],[preselected]);
%
%   this function is deprecated, use analysesList = selectInList(thisView,'analyses',title,preselected)
%
% $Id$ 

analysesList = selectInList(thisView,'analyses',varargin);
mrWarnDlg('(selectAnalyses) selectAnalyses is deprecated. Please use ''selectInList(view,''analyses'',...)'' instead');