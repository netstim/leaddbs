% getAlignGUI.m
%
%        $Id$ 
%      usage: [trans rot] = getAlignGUI(handles)
%         by: justin gardner
%       date: 05/16/09
%    purpose: returns the trans and rot of what the GUI is set to. Used for left/right arrow buttons
%
function [trans rot] = getAlignGUI(handles)

% check arguments
if ~any(nargin == [1])
  help getAlignGUI
  return
end

trans = [str2double(get(handles.transX,'String')),...
    str2double(get(handles.transY,'String')),...
    str2double(get(handles.transZ,'String'))];
% str2double returns nan for non-numeric strings, change to zeros
nanIndices = find(isnan(trans));
trans(nanIndices) = zeros(size(nanIndices));
		
% Get rot from GUI
rot = [str2double(get(handles.rotX,'String')),...
    str2double(get(handles.rotY,'String')),...
    str2double(get(handles.rotZ,'String'))];
% str2double returns nan for non-numeric strings, change to zeros
nanIndices = find(isnan(rot));
rot(nanIndices) = zeros(size(nanIndices));

