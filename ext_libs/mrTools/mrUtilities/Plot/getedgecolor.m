% getedgecolor.m
%
%      usage: getedgecolor()
%         by: justin gardner
%       date: 03/02/06
%    purpose: This is called to get the edge color for
%             graphs. e.g. for exporting figures we 
%             sometimes want edges to be white rather
%             than black and this is how functions
%             know what color to plot in.
%
function [edgecolor othercolors] = getedgecolor()

% check arguments
if ~any(nargin == [0])
  help getedgecolor
  return
end

% get the edgecolor, this should be
% either 'k' or white = [0.99 0.99 0.99]
% if it is set to 'w' it will be reversed
% by matlab to black when you export the figure
global edgecolor;
if (isempty(edgecolor))
  edgecolor = 'k';
end
if (isequal(edgecolor,[0.99 0.99 0.99]))
  set(gcf,'Color','k');
  othercolors = 'gwrycm';
else
  othercolors = 'gkrbcm';
end


