% mrMoveWindow.m
%
%      usage: mrMoveWindow(<newpos>)
%         by: justin gardner
%       date: 08/08/07
%    purpose: moves mrLoadRet window. This is useful if the
%             mrLoadRet window gets auto placed somewhere 
%             where it is not visible (like on a second monitor
%             that is no longer visible). If newpos is provided
%             it can be array of length 2 specifying [x y]
%             or an array of length 4 specifying [x y w h]. 
%             Default is to move the window to the bottom left
%             corner of the primary display without changing
%             with and height--if it is deemed to be off the
%             primary display.
% 
%
function retval = mrMoveWindow(newpos)

% check arguments
if ~any(nargin == [0 1])
  help mrMoveWindow
  return
end

% min and max positions that should be on primary screen
minXpos = 0;
maxXpos = 1200;
minYpos = 0;
maxYpos = 800;

% default positions to change to
defaultXpos = 1;
defaultYpos = 1;

% get mrGlobals
mrGlobals

for viewNum = 1:length(MLR.views)
  if isview(MLR.views{viewNum})
    % get this view
    v = MLR.views{viewNum};
    if ishandle(viewGet(v,'fignum'))
      % get its screen position and print
      pos = get(viewGet(v,'fignum'),'Position');
      disp(sprintf('(mrMoveWindow) View %i is at [%i %i %i %i]',viewNum,pos(1),pos(2),pos(3),pos(4)));
      % if it is out of bounds, then move it
      if ieNotDefined('newpos')
	outOfBounds = 0;
	if ((pos(1)+pos(3)) < minXpos) || (pos(1) > maxXpos)
	  outOfBounds = 1;
	  pos(1) = defaultXpos;
	end
	if ((pos(2)+pos(4)) < minYpos) || (pos(2) > maxYpos)
	  outOfBounds = 1;
	  pos(2) = defaultYpos;
	end
      else
	% set the position according to passed in values
	outOfBounds = 1;
	if length(newpos) == 2
	  pos(1:2) = newpos(1:2);
	elseif length(newpos) == 4
	  pos = newpos;
	else
	  disp(sprintf('(mrMoveWindow) Newpos not length of 2 or 4'));
	end
      end
      % now actually move
      if outOfBounds
	disp(sprintf('(mrMoveWindow) View %i has been moved to [%i %i %i %i]',viewNum,pos(1),pos(2),pos(3),pos(4)));
	set(viewGet(v,'fignum'),'Position',pos);
	refreshMLRDisplay(viewNum);
      end
    end
  end
end
