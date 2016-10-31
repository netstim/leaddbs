% makeEqualXaxis.m
%
%      usage: makeEqualXaxis(nrows,ncols,<axisNums>)
%         by: justin gardner
%       date: 05/14/08
%    purpose: Make equal x axis across all subplots.
%             Pass in the number of rows, columns of
%             subplots and then all the axis numbers that
%             you want to equalize over.
%
function makeEqualXaxis(nrows,ncols,axisNums)

% check arguments
if ~any(nargin == [2 3])
  help makeEqualXaxis
  return
end

% default axisNums
if ieNotDefined('axisNums')
  axisNums = 1:(nrows*ncols);
end

% get the smallest and largest yaxis
xmin = inf;xmax = -inf;
for aNum = axisNums
  subplot(nrows,ncols,aNum);
  a = axis;
  xmin = min(xmin,a(1));
  xmax = max(xmax,a(2));
end

% and set them
for aNum = axisNums
  subplot(nrows,ncols,aNum);
  xaxis(xmin,xmax);
end

