% getMonitorNumber.m
%
%      usage: monitorNumbers = getMonitorNumber(figurePositions,monitorPositions)
%         by: julien besle
%       date: 19/12/2010
%        $Id$
%    purpose: finds which monitor most of a figure (or several figures) is diplayed on, 
%             or the closest monitor if the figure is outside any monitor area
%             optionally returns a new position Figure position

function varargout = getMonitorNumber(figurePositions,monitorPositions)

if nargout==2
  varargout = cell(1,2);
else
  varargout{1} = [];
end

minIntersection = 150; %minimum intersection of the figure with the identified monitor for new figure position(in pixels)

%convert to [left bottom];[right top]
figureLeftBottom = figurePositions(:,1:2);
figureRigthTop = figurePositions(:,3:4)+figurePositions(:,1:2);
monitorLeftBottom = monitorPositions(:,1:2);
monitorRigthTop = monitorPositions(:,3:4)+monitorPositions(:,1:2);

nFigures = size(figurePositions,1);
nMonitors = size(monitorPositions,1);
monitorNumbers = zeros(1,nFigures);
if nargout==2
  newFigurePositions = figurePositions;
end

for iFigure = 1:size(figurePositions,1)
  intersectionDims =  min(repmat(figureRigthTop(iFigure,:),nMonitors,1),monitorRigthTop) - ...
    max(repmat(figureLeftBottom(iFigure,:),nMonitors,1),monitorLeftBottom);
  if any(all(intersectionDims>0,2)) %if any intersection is positive on both horizontal and vertical axes
    %the monitor is the one with the largest intersection area
    thisIntersectionDims = intersectionDims; 
    thisIntersectionDims(intersectionDims<0) = 0; %Negative intersections are considered 0 
    [maxIntersection,monitorNumbers(iFigure)] = max(prod(thisIntersectionDims,2));
  else %otherwise we have to find which monitor is closest
    %we'll just take the intersection with the smallest distance on any axis (vertical,horizontal or diagonal)
    %that is the length of the diagonal for double negative intersections
    %or, for intersections that are negative only on one axis, the distance on this axis 
    thisIntersectionDims = intersectionDims; 
    thisIntersectionDims(intersectionDims>0) = 0; %Positive distances are considered 0 (any positive distance is an overlap in a given axis)
    [minDistance,monitorNumbers(iFigure)] = min(sqrt(sum(thisIntersectionDims.^2,2)));
    
  end
  %if asked for, compute new positions for figure
  if nargout==2
    %correct left position if intersection negative on horizontal axis
    if intersectionDims(monitorNumbers(iFigure),1)<minIntersection
      if figurePositions(iFigure,1)<monitorPositions(monitorNumbers(iFigure),1)
        newFigurePositions(iFigure,1) = figurePositions(iFigure,1)-intersectionDims(monitorNumbers(iFigure),1)+... %shift figure right to edge of monitor
          minIntersection; %add a hundred pixels
%             min(figurePositions(iFigure,3),monitorPositions(monitorNumbers(iFigure),3)); %and shift more to see everything that fits on the monitor
      else
        newFigurePositions(iFigure,1) = figurePositions(iFigure,1)+intersectionDims(monitorNumbers(iFigure),1)-... %shift figure left to edge of monitor
          minIntersection; %remove a hundred pixels
      end
    end
    %correct bottom position if intersection negative on vertical axis
    if intersectionDims(monitorNumbers(iFigure),2)<minIntersection
      if figurePositions(iFigure,2)<monitorPositions(monitorNumbers(iFigure),2)
        newFigurePositions(iFigure,2) = figurePositions(iFigure,2)-intersectionDims(monitorNumbers(iFigure),2)+... %shift figure up to edge of monitor
          minIntersection; %add a hundred pixels
      else
        %if the figure is too high, we want to make sure we're able to grab the top bar
        newFigurePositions(iFigure,2) = figurePositions(iFigure,2)+intersectionDims(monitorNumbers(iFigure),2)-... %shift figure down to edge of monitor
           figurePositions(iFigure,4)-30; %and shift more to make sure we can see the top bar
      end
    end
  end
end
 
varargout{1} = monitorNumbers;
if nargout==2
  varargout{2} = newFigurePositions;
end


    