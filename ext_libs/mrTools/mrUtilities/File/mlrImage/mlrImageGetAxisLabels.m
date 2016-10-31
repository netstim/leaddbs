% mlrImageGetAxisLabels.m
%
%        $Id:$ 
%      usage: axisLabels = mlrImageGetAxisLabels(xform)
%         by: justin gardner
%       date: 08/19/11
%    purpose: This function uses the passed in xform (e.g. a qform) to determine in which direction
%             each axis of the image goes. axisLabels contains the follwing fields:
%
%             labels: are readable labels
%             dirLabels: is a cell array with two strings for each direction (-/+) for each axis.
%             mapping: is the map of each axis to the closest axis in the magnet space. That is
%               each element is the axis in the image space and what axis it corresponds to in the
%               magnet space.
%             reverseMapping: is the opposite of mapping
%             dir: is the direction in which that axis goes in the magnet space.
%             orient: is a 3 letter string representing the axis orientation (like LPI)
%
function retval = mlrImageGetAxisLabels(xform)

% default return argument
retval = [];

% check input arguments
if nargin ~= 1
  help mlrImageGetAxisLabels;
  return
end

% check for empty xform
if isempty(xform)
  disp(sprintf('(mlrImageGetAxisLabels) Empty xform'));
  return
end

axisNames  = {'X','Y','Z'};
for axisNum = 1:3
  % get the vector of the axis in the image that we want to label
  axisVector = zeros(3,1);
  axisVector(axisNum) = 1;

  % find out which direction in magnet space the axis vector goes
  axisDirVector = xform(1:3,1:3)*axisVector;

  % normalize to unit length
  axisDirVector = axisDirVector ./ sqrt(sum(axisDirVector.^2));

  % these are the magnet cardinal axis and their names
  cardinalAxisDirs = {[1 0 0],[0 1 0],[0 0 1],[-1 0 0],[0 -1 0],[0 0 -1]};
  cardinalAxisLabels = {'right','anterior','superior','left','posterior','inferior'};

  % now get the angle of this axisDirVector with
  % each of the magnet cardinal axis. Note that
  % we are computing for both directions of the axis
  % for conveinence (i.e. left and right)
  for i = 1:length(cardinalAxisDirs)
    angles(i) = r2d(acos(dot(axisDirVector,cardinalAxisDirs{i})));
  end

  % sort the angles (remembering which axis they originally came from). Thus
  % sortedAxisNum contains an ordered list of which axis the vector is closest
  % too. The closest axis is sortedAxisNum(1) and the farthest axis is sortedAxisNum(6)
  [angles sortedAxisNum] = sort(angles);

  % get the closest axis direction (i.e. the one with the smallest angle which
  % is the first in the list of sortedAxisNum)
  axisDirs(axisNum,:) = cardinalAxisDirs{sortedAxisNum(1)}';
  
  % if the closest angle is less than an arbitrary value than we will consider
  % the axis to be a pure direction - if not, we will label
  % with a combination of the two closest axis.
  if angles(1) < 5
    axisDirLabels{axisNum} = {cardinalAxisLabels{sortedAxisNum(6)} cardinalAxisLabels{sortedAxisNum(1)}};
  else
    axisDirLabels{axisNum} = {sprintf('%s/%s',cardinalAxisLabels{sortedAxisNum(6)},cardinalAxisLabels{sortedAxisNum(5)}) sprintf('%s/%s',cardinalAxisLabels{sortedAxisNum(1)},cardinalAxisLabels{sortedAxisNum(2)})};
  end
  % make a single name for each axis
  axisLabels{axisNum} = sprintf('%s <- %s -> %s',axisDirLabels{axisNum}{1},axisNames{axisNum},axisDirLabels{axisNum}{2});
end

% convert the axisDirs to axisMapping (i.e. what axis in the magnet each axis in the image
% corresponds to and in which direction it points
[axisReverseMapping row axisDir] = find(axisDirs);
[dummy axisMapping] = sort(axisReverseMapping);
axisDir = axisDir(axisMapping);

% get the orientation string (like LPI)
orient = sprintf('%s%s%s',upper(axisDirLabels{1}{1}(1)),upper(axisDirLabels{2}{1}(1)),upper(axisDirLabels{3}{1}(1)));

% return structure with all ifno
retval.labels = axisLabels;
retval.dirLabels = axisDirLabels;
retval.mapping = axisMapping;
retval.reverseMapping = axisReverseMapping;
retval.dir = axisDir;
retval.orient = orient;

% convert radians to degrees
%
% usage: degrees = r2d(radians);
function degrees = r2d(angle)

degrees = (angle/(2*pi))*360;

% if larger than 360 degrees then subtract
% 360 degrees
while (sum(degrees>360))
  degrees = degrees - (degrees>360)*360;
end

% if less than 360 degreees then add 
% 360 degrees
while (sum(degrees<-360))
  degrees = degrees + (degrees<-360)*360;
end
