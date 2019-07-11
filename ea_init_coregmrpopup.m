function ea_init_coregmrpopup(handles,refine)
% 
%
% USAGE:
%
%    ea_init_coregmrpopup(handles,refine)
%
% INPUTS:
%    handles:
%    refine:
%
% .. AUTHOR:
%       - Andreas Horn, Original file
%       - Ningfei Li, Original file
%       - Daniel Duarte, Documentation

if ~exist('refine','var')
    refine=0;
end

cmethods={'SPM',...
    'FSL FLIRT',...
    'ANTs',...
    'BRAINSFIT',...
    'Hybrid SPM & ANTs',...
    'Hybrid SPM & FSL',...
    'Hybrid SPM & BRAINSFIT'};

set(handles.coregmrpopup,'String',cmethods)
set(handles.coregmrpopup,'Value',1); % default SPM
