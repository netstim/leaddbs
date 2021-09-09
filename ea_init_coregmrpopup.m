function ea_init_coregmrpopup(handles)
% Initializes MR coregistration methods popupmenu.
%
% USAGE:
%
%    ea_init_coregmrpopup(handles,refine)
%
% INPUTS:
%    handles:       LEAD GUI handle
%    refine:        deprecated

cmethods={'SPM',...
    'FSL FLIRT',...
    'ANTs',...
    'BRAINSFIT',...
    'Hybrid SPM & ANTs',...
    'Hybrid SPM & FSL',...
    'Hybrid SPM & BRAINSFIT'};

set(handles.coregmrmethod, 'String', cmethods);
set(handles.coregmrmethod, 'Value', 1); % default SPM
