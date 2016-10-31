% correctMonitorPosition.m
%
%      usage: monitorPositions = correctMonitorPosition(monitorPositions)
%         by: julien besle
%       date: 19/12/2010
%        $Id$
%    purpose: correct MonitorPosition property or root object
%             returned by get(0,'monitorPositions')
%             so that additional monitor positions are in the form [left bottom width height]
%               with the horizontal axis increasing from left to right
%               and the vertical axis increasing from bottom to top
%               with an origin at the left/bottm of the primary monitor
%             consistent with the usual position vector form 
%            (as used by the 'position' property of figures)

function monitorPositions = correctMonitorPosition(monitorPositions)

if size(monitorPositions,1)>1
  if ispc 
    %windows machines return the position of secondary monitors in a 
    % [left top right bottom] form, where
    %       - the vertical origin is the top of the PRIMARY monitor
    %       - the horizontal origin is the left of the PRIMARY monitor
    %       - the vertical axis increases from top to bottom
    %       - the origin of the PRIMARY monitor (the top left) is (1,1) 
    %         (but we don't care about this because it's only one pixel)

  %the first position vector is the primary monitor in the correct forms. don't touch it
  % change the origin of the vertical axis, flip it and swap top and bottom
  monitorPositions(2:end,[4 2]) = (monitorPositions(2:end,[2 4]) - monitorPositions(1,4))*-1 +1;

  elseif ismac
    %macs return the position of all monitors in a 
    % [left top width heigth] form, where
    %       - the vertical origin is the top of the TOP monitor
    %       - the horizontal origin is the left of the LEFTMOST monitor
    %       - the vertical axis increases from top to bottom
    %       - the origin of the LEFTMOST,TOP monitor (the top left) is (1,0) 
    %         (but we don't care about this because it's only one pixel)

  %first we need to change the position to an origin at the left and bottom of the PRIMARY monitor (the coordinate system of figure positions)
  % put origin to the left of the PRIMARY monitor (first line)
  monitorPositions(:,1) = monitorPositions(:,1) - monitorPositions(1,1);
  % shift the vertical origin to the bottom of the topmost monitor
  monitorPositions(:,2) = monitorPositions(:,2) + monitorPositions(monitorPositions(:,2)==0,4);
  %put origin to the bottom of the PRIMARY monitor
  monitorPositions(:,2) = monitorPositions(:,2) - monitorPositions(1,2);
 
    
  else
    keyboard
    %I don't know about other platforms... if they do what's described in the Matlab documentation
    %we shouldn't have to do anything
  end
  
  % compute width and heigth
  monitorPositions(2:end,3:4) = monitorPositions(2:end,3:4)+1 - monitorPositions(2:end,1:2);
end