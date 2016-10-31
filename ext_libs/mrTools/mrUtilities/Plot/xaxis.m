
%
%      usage: xaxis(xmin,xmax)
%         by: justin gardner
%       date: 02/17/03
%    purpose: sets the x axis of a figure. With no arguments returns xmin and xmax that is currently set
%
function [xmin xmax] = xaxis(xmin,xmax)

if ((nargin == 1) && (length(xmin) == 2))
  xmax = xmin(2);
  xmin = xmin(1);
elseif (nargin == 0)
  a = axis;
  xmin = a(1);
  xmax = a(2);
  return
elseif (nargin ~= 2)
  help xaxis;
  return
end

a = axis;

% default to what is already there
if isempty(xmin),xmin = a(1);end
if isempty(xmax),xmax = a(2);end

if xmin<xmax
  axis([xmin xmax a(3) a(4)]);
end
