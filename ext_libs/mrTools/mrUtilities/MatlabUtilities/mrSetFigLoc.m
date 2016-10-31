function mrSetFigLoc(figname,pos)
%
% mrSetFigLoc(figname,[pos])
%
% Sets a field in the global variable mrDEFAULTS.figloc, which is a
% structure with fields for each figure name.
%
% figname is a string that specifies each type of figure window
% pos is a 4-vector specifying lowerleft corner and size
%     default: [100 100 560 420]
%
% Examples:
%   mrSetFigLoc('mrLoadRetGUI',[100 100 560 420]);
%
% djh, 5/2007

global mrDEFAULTS

% if mrDEFAULTS is empty, then we should try to load it
if isempty(mrDEFAULTS)
  mrDEFAULTS = loadMrDefaults;
end

if ieNotDefined('pos')
    pos = [100 100 560 420];
end

% mrDEFAULTS.figloc = setfield(mrDEFAULTS.figloc,figname,pos); 
% crashes on MATLAB Version 7.4.0.287 (R2007a) if figloc field does not exist

mrDEFAULTS.figloc.(figname) = pos; % dynamic field assignment creats field if not present
saveMrDefaults