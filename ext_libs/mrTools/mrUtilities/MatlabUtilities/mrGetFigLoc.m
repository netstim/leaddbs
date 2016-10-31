function pos = mrGetFigLoc(figname)
%
% pos = mrGetFigLoc(figname)
%
% Gets a field in the global variable mrDEFAULTS.figloc, which is a
% structure with fields for each figure name.
%
% figname is a string that specifies each type of figure window
% pos is a 4-vector specifying lowerleft corner and size
%
% Examples:
%   pos = mrGetFigLoc('mrLoadRetGUI');
%   pos = mrGetFigLoc('buttondlg');
%   pos = mrGetFigLoc('graphFigure');
%   pos = mrGetFigLoc('mrParamsDialog');
%   pos = mrGetFigLoc('mrParamsDialogHelp');
%
% djh, 5/2007

global mrDEFAULTS

% if mrDEFAULTS is empty, then we should try to load it
if isempty(mrDEFAULTS)
  mrDEFAULTS = loadMrDefaults;
end

if ~isempty(mrDEFAULTS) && isfield(mrDEFAULTS.figloc,figname)
    % pos = getfield(mrDEFAULTS.figloc,figname);
    pos = mrDEFAULTS.figloc.(figname);
else
    pos = [];
end

if isunix && strcmp(version,'7.14.0.739 (R2012a)') %some bugs in that version, disappears in next version, 
                                                   %don't know if same happens in previous versions
  pause(.1);                          % basically there is a mismatch between outerposition and position that 
                                      % messes up the coordinates inside the figure but disappears when
                                      % entering debug mode or pausing
%   figloc = figloc + [0 38 -16 -38];  % if not pausing, this is the difference in location/size values   
%   disp(['Position ' mat2str(get(h,'position'))]);
end

if ~verLessThan('matlab','8.4') && ~ispc && ~isunix && ~isempty(pos)  
  pos(2) = pos(2)+30;
end
  
