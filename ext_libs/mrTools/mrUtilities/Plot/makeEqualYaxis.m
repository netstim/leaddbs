% makeEqualYaxis.m
%
%      usage: makeEqualYaxis(nrows,ncols,<axisNums>)
%         by: justin gardner
%       date: 05/14/08
%    purpose: Make equal y axis across all subplots.
%             Pass in the number of rows, columns of
%             subplots and then all the axis numbers that
%             you want to equalize over.
%
function makeEqualYaxis(nrows,ncols,axisNums)

% check arguments
if ~any(nargin == [2 3])
  help makeEqualYaxis
  return
end

% default axisNums
if ieNotDefined('axisNums')
  axisNums = 1:(nrows*ncols);
end

% get the smallest and largest yaxis
ymin = inf;ymax = -inf;
for aNum = axisNums
  subplot(nrows,ncols,aNum);
  a = axis;
  ymin = min(ymin,a(3));
  ymax = max(ymax,a(4));
end

% and set them
for aNum = axisNums
  subplot(nrows,ncols,aNum);
  yaxis(ymin,ymax);
end

