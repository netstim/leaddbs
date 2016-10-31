% yaxis
%
%      usage: yaxis(ymin,ymax)
%         by: justin gardner
%       date: 02/17/03
%    purpose: sets the y axis of a figure. With no arguments returns xmin and xmax that is currently set
%
function [ymin ymax] = yaxis(ymin,ymax)

if ((nargin == 1) && (length(ymin) == 2))
  ymax = ymin(2);
  ymin = ymin(1);
elseif (nargin == 0)
  a = axis;
  ymin = a(3);
  ymax = a(4);
  return
elseif (nargin ~= 2)
  help yaxis;
  return
end

a = axis;

% default to what is already there
if isempty(ymin),ymin=a(3);end
if isempty(ymax),ymax=a(4);end

if ymin < ymax
  axis([a(1) a(2) ymin ymax]);
end
