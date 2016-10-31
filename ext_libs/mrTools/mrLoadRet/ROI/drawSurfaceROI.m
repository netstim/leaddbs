% drawSurfaceROI.m
%
%        $Id$ 
%      usage: coords = drawSurfaceROI(v)
%         by: Luke Robot Hospadaruk <hospada1@msu.edu>
%       date: 11/28/09
%    purpose: draw an ROI on 3D surface
%
function coords = drawSurfaceROI(v)

coords = [];
%use modified getpts function to grab the correct screen locations

% get figure and make sure hold is on
fig = viewGet(v,'figureNumber');
hold(get(fig, 'CurrentAxes'), 'on');
roiKeyVertices = [];
% for surfaces in matlab past 8.4 select3d is busted, so need a new
% method to get points which is implemented mlrGetptsSurface
if ~verLessThan('matlab','8.4') && (viewGet(v,'baseType') == 2)
  % tell the user what to do
  oneTimeWarning('drawROIOnSurface','(drawROI) Click on each boundary vertex of ROI. To finish making the ROI control-click (click while holding down the control key) on a point that you want inside the ROI boundary',-1);
  % get the display surface
  h = viewGet(v,'baseHandle');
  if ~isempty(h)
    % and get user points that are clicked on
    roiKeyVertices = mlrGetptsSurface(v,h);
  else
    disp(sprintf('(drawSurfaceROI) Could not find surface'));
    return
  end
else
  % tell the user what to do
  oneTimeWarning('drawROIOnSurface','(drawROI) Click on each boundary vertex of ROI. To finish making the ROI double click on a point that you want inside the ROI boundary',-1);
  % old style (before mathworks broke select3d)
  [xi, yi] = mlr_getpts_3d(fig);
  % calculate key vertices from list of x, y points
  disp('(drawSurfaceROI) Calculating key vertices');
  for i=1:(numel(xi))
    [x y s xBase yBase sBase xT yT sT vi] = getMouseCoords(viewGet(v,'viewNum'), [xi(i) yi(i)]);
    roiKeyVertices = [roiKeyVertices vi];
  end
end

%make sure there are at least 3 perimiter and 1 interior vertices
if(numel(roiKeyVertices) < 4)
  mrWarnDlg('(drawSurfaceROI) Not enough vertices selected, you need at least 4');
  coords = [];
  return
end

% get coord map for baseqw
cmap = viewGet(v,'baseCoordMap');

% get all vertices in boundary
disp('(drawSurfaceROI) Setting up boundary search');
% get the base surface
baseSurface = viewGet(v,'baseSurface');
% set up the list of vertices to go around back to the first in the list
drawVertices = [roiKeyVertices(1:end-1) roiKeyVertices(1)];
% list of all edge vertices
edgeVertices = [];
for vertexIndex = 2:length(drawVertices)
  pathVertices = mlrGetPathBetween(baseSurface,drawVertices(vertexIndex-1),drawVertices(vertexIndex));
  % keep list of all vertices
  edgeVertices = [edgeVertices pathVertices];
end
% plot the edge
plot3(baseSurface.vtcs(edgeVertices,1),baseSurface.vtcs(edgeVertices,2),baseSurface.vtcs(edgeVertices,3),'r.');

% start list
vList = roiKeyVertices(end);
% get surface
baseSurface = viewGet(v,'baseSurface');
tris = baseSurface.tris(:);
trisSize = size(baseSurface.tris);
lastLen = -inf;
% find all neighboring vertices, by getting the location
% of all current vertex in (linearized) tris list
while(length(vList)>lastLen)
  % keep the last length to make sure we are growing
  lastLen = length(vList);
  % find all the triangles that intersect with the current list
  intersectTris = find(ismember(tris,vList));
  [intersectTris,dummy] = ind2sub(trisSize,intersectTris);
  intersectTris = unique(intersectTris);
  % now we know all triangles that intersect so find all neighboring points
  oldVList = vList;
  vList = unique(tris(sub2ind(trisSize,repmat(intersectTris,1,3),repmat([1 2 3],length(intersectTris),1))));
  % remove edge vertices
  vList = setdiff(vList,edgeVertices);
  % plot
  drawList = setdiff(vList,oldVList);
  plot3(baseSurface.vtcs(drawList,1),baseSurface.vtcs(drawList,2),baseSurface.vtcs(drawList,3),'w.');
  drawnow;
end
% add back edge vertices
vList = union(vList,edgeVertices);
hold(get(fig, 'CurrentAxes'), 'off');

coords = [];
% cycle over cortical depths and add points
corticalDepths = viewGet(v,'corticalDepth');
for corticalDepth = corticalDepths(1):0.1:corticalDepths(end)
  % get the cortical slice index
  corticalSlice = find(cmap.corticalDepths == (round(corticalDepth*10)/10));
  % read off coordinates for this depth
  coords = [coords;squeeze(cmap.coords(1,vList,1,:,corticalSlice))];
end
% return unique coords
coords = unique(coords,'rows');
coords = [coords ones(size(coords,1),1)]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get current mouse position in image coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x y s xBase yBase sBase xTal yTal zTal vi] = getMouseCoords(viewNum, point)

% get the view
mrGlobals;
view = viewGet([],'view',viewNum);

% no bases
if viewGet(view,'numBase') == 0
	x = nan;y= nan;s = nan;xBase = nan;yBase = nan;sBase = nan;xTal = nan;yTal = nan;zTal = nan;
	return
end
baseType = viewGet(view,'baseType');


% get location of pointer
pointerLoc = point;
%set it as the current point so select3d will get the right point
set(view.figure, 'CurrentPoint', point);
%pointerLoc = get(MLR.interrogator{viewNum}.axesnum,'CurrentPoint');

if baseType <= 1
	mouseY = round(pointerLoc(1,1));
	mouseX = round(pointerLoc(1,2));
	
	% get base coordinates
	baseCoords = viewGet(view,'cursliceBaseCoords');
	% convert mouse to baseCoords
	if (mouseX>0) && (mouseX<=size(baseCoords,1)) && (mouseY>0) && (mouseY<=size(baseCoords,2))
		xBase = baseCoords(mouseX,mouseY,1);
		yBase = baseCoords(mouseX,mouseY,2);
		sBase = baseCoords(mouseX,mouseY,3);
	else
		x = nan;y = nan; s = nan;
		xBase = nan;yBase = nan; sBase = nan;
		xTal = nan; yTal = nan; zTal = nan;
		return
	end
else
	% handle getting coordinates for surface
	baseSurface = viewGet(view,'baseSurface');
	baseDims = viewGet(view,'baseSurfaceDims');
	pos = [];xBase = nan; yBase = nan; sBase = nan;
	% check mouse bounding box coords against baseDims
	% for a quick check to see if we are in the volume
	if all(pointerLoc(1,:) >= 0)
		% then use select3d which is slooow, but accurate
		hobj = get(get(view.figure, 'CurrentAxes'),'Children');
		% make sure we are using the correct object (should be the 3D
		% brain). Use end here because with searchForVoxel we plot a
		% point on the image and the brain object always seems to be
		% last. But if this is not always the case, then we may need
		% to do a little more work here to find the correct object
		[pos v vi] = select3d(hobj(end));
		% convert the index to the coordinates
		if ~isempty(pos)
			baseCoordMap = viewGet(view,'baseCoordMap');
      %we'll take the coordinates of the middle of whatever range of cortical depth is currenlty selected
      corticalSlice = ceil(mean(viewGet(view,'corticalDepth'))*size(baseCoordMap.coords,5));
			pos = round(squeeze(baseCoordMap.coords(1,vi,1,:,corticalSlice)));
			xBase = pos(1);yBase = pos(2);sBase = pos(3);
		end
	end
end

% transform from base coordinates to talairach coordinates, if the base has a talairach transform defined
base2tal = viewGet(view,'base2tal'); % keyboard
if(~isempty(base2tal))
	talCoords = round(base2tal * [xBase yBase sBase 1]');
	xTal = talCoords(1); yTal = talCoords(2); zTal = talCoords(3);
else
	xTal = nan; yTal = nan; zTal = nan;
end


% transform from base coordinates into scan coordinates
base2scan = viewGet(view,'base2scan');
if isempty(base2scan), x = nan; y=nan; s=nan; return,  end
transformed = round(base2scan*[xBase yBase sBase 1]');

x = transformed(1);
y = transformed(2);
s = transformed(3);

% get the scan dims to make sure we haven't jumped off end
% of scan
scanDims = viewGet(view,'scanDims',viewGet(view,'curScan'));
if ((x < 1) || (x > scanDims(1)) || ...
		(y < 1) || (y > scanDims(2)) || ...
		(s < 1) || (s > scanDims(3)))
	x = nan;y = nan;s = nan;
end


function [x,y] = mlr_getpts_3d(varargin)
%Modified version to get points for getMouseCoords from 3d surface
%for ROI creation on surfaces
%GETPTS Select points with mouse.
%   [X,Y] = GETPTS(FIG) lets you choose a set of points in the
%   current axes of figure FIG using the mouse. Coordinates of
%   the selected points are returned in the vectors X and Y. Use
%   normal button clicks to add points.  A shift-, right-, or
%   double-click adds a final point and ends the selection.
%   Pressing RETURN or ENTER ends the selection without adding
%   a final point.  Pressing BACKSPACE or DELETE removes the
%   previously selected point.
%
%   [X,Y] = GETPTS(AX) lets you choose points in the axes
%   specified by the handle AX.
%
%   [X,Y] = GETPTS is the same as [X,Y] = GETPTS(GCF).
%
%   Example
%   --------
%       imshow('moon.tif')
%       [x,y] = getpts_luke
%
%   See also GETRECT, GETLINE.

%   Callback syntaxes:
%       getpts_luke('KeyPress')
%       getpts_luke('FirstButtonDown')
%       getpts_luke('NextButtonDown')

%   Copyright 1993-2008 The MathWorks, Inc.
%   $Revision: 1533 $  $Date: 2009-10-07 02:26:45 +0900 (Wed, 07 Oct 2009) $
%	Modified by: Luke Hospadaruk

global GETPTS_FIG GETPTS_AX GETPTS_H1 GETPTS_H2
global GETPTS_PT1 line_points
if ((nargin >= 1) && (ischar(varargin{1})))
	% Callback invocation: 'KeyPress', 'FirstButtonDown', or
	% 'NextButtonDown'.
	feval(varargin{:});
	return;
end

if (nargin < 1)
	GETPTS_AX = gca;
	GETPTS_FIG = ancestor(GETPTS_AX, 'figure');
else
	if (~ishandle(varargin{1}))
		eid = 'Images:getpts:expectedHandle';
		error(eid, '%s', 'First argument is not a valid handle');
	end
	
	switch get(varargin{1}, 'Type')
		case 'figure'
			GETPTS_FIG = varargin{1};
			GETPTS_AX = get(GETPTS_FIG, 'CurrentAxes');
			if (isempty(GETPTS_AX))
				GETPTS_AX = axes('Parent', GETPTS_FIG);
			end
			
		case 'axes'
			GETPTS_AX = varargin{1};
			GETPTS_FIG = ancestor(GETPTS_AX, 'figure');
			
		otherwise
			eid = 'Images:getpts:expectedFigureOrAxesHandle';
			error(eid, '%s', 'First argument should be a figure or axes handle');
			
	end
end

% Bring target figure forward
figure(GETPTS_FIG);

% Remember initial figure state
state = uisuspend(GETPTS_FIG);

% Set up initial callbacks for initial stage
[pointerShape, pointerHotSpot] = CreatePointer;
set(GETPTS_FIG, 'WindowButtonDownFcn', 'mlr_getpts_3d(''FirstButtonDown'');', ...
	'KeyPressFcn', 'mlr_getpts_3d(''KeyPress'');', ...
	'Pointer', 'custom', ...
	'PointerShapeCData', pointerShape, ...
	'PointerShapeHotSpot', pointerHotSpot);

% Initialize the lines to be used for the drag
markerSize = 9;
GETPTS_H1 = line('Parent', GETPTS_AX, ...
	'XData', [], ...
	'YData', [], ...
	'Visible', 'off', ...
	'Clipping', 'off', ...
	'Color', 'c', ...
	'LineStyle', 'none', ...
	'Marker', '+', ...
	'MarkerSize', markerSize);

GETPTS_H2 = line('Parent', GETPTS_AX, ...
	'XData', [], ...
	'YData', [], ...
	'Visible', 'off', ...
	'Clipping', 'off', ...
	'Color', 'm', ...
	'LineStyle', 'none', ...
	'Marker', 'x', ...
	'MarkerSize', markerSize);

% We're ready; wait for the user to do the drag
% Wrap the call to waitfor in try-catch so we'll
% have a chance to clean up after ourselves.
errCatch = 0;
try
	waitfor(GETPTS_H1, 'UserData', 'Completed');
catch
	errCatch=1;
end

% After the waitfor, if GETPTS_H1 is still valid
% and its UserData is 'Completed', then the user
% completed the drag.  If not, the user interrupted
% the action somehow, perhaps by a Ctrl-C in the
% command window or by closing the figure.

if (errCatch == 1)
	errStatus = 'trap';
	
elseif (~ishandle(GETPTS_H1) || ...
		~strcmp(get(GETPTS_H1, 'UserData'), 'Completed'))
	errStatus = 'unknown';
	
else
	errStatus = 'ok';
	x = get(GETPTS_H1, 'XData');
	y = get(GETPTS_H1, 'YData');
	x = x(:);
	y = y(:);
	% If no points were selected, return rectangular empties.
	% This makes it easier to handle degenerate cases in
	% functions that call getpts_luke.
	if (isempty(x))
		x = zeros(0,1);
	end
	if (isempty(y))
		y = zeros(0,1);
	end
end

% Delete the animation objects
if (ishandle(GETPTS_H1))
	delete(GETPTS_H1);
end
if (ishandle(GETPTS_H2))
	delete(GETPTS_H2);
end

% Restore the figure state
if (ishandle(GETPTS_FIG))
	uirestore(state);
end

% Clean up the global workspace
clear global GETPTS_FIG GETPTS_AX GETPTS_H1 GETPTS_H2
clear global GETPTS_PT1 line_points

% Depending on the error status, return the answer or generate
% an error message.
switch errStatus
	case 'ok'
		% No action needed.
		
	case 'trap'
		% An error was trapped during the waitfor
		eid = 'Images:getpts:interruptedMouseSelection';
		error(eid, '%s', 'Interruption during mouse point selection.');
		
	case 'unknown'
		% User did something to cause the point selection to
		% terminate abnormally.  For example, we would get here
		% if the user closed the figure in the middle of the selection.
		eid = 'Images:getpts:interruptedMouseSelection';
		error(eid, '%s', 'Interruption during mouse point selection.');
end


%--------------------------------------------------
% Subfunction KeyPress
%--------------------------------------------------
function KeyPress %#ok

global GETPTS_FIG GETPTS_AX GETPTS_H1 GETPTS_H2
global GETPTS_PT1 line_points

key = get(GETPTS_FIG, 'CurrentCharacter');
switch key
	case {char(8), char(127)}  % delete and backspace keys
		x = get(GETPTS_H1, 'XData');
		y = get(GETPTS_H1, 'YData');
		switch length(x)
			case 0
				% nothing to do
			case 1
				% remove point and start over
				set([GETPTS_H1 GETPTS_H2], ...
					'XData', [], ...
					'YData', []);
				set(GETPTS_FIG, 'WindowButtonDownFcn', ...
					'getpts_luke(''FirstButtonDown'');');
			otherwise
				% remove last point
				set([GETPTS_H1 GETPTS_H2], ...
					'XData', x(1:end-1), ...
					'YData', y(1:end-1));
		end
		
	case {char(13), char(3)}   % enter and return keys
		% return control to line after waitfor
		set(GETPTS_H1, 'UserData', 'Completed');
		
end

%--------------------------------------------------
% Subfunction FirstButtonDown
%--------------------------------------------------
function FirstButtonDown %#ok

global GETPTS_FIG GETPTS_AX GETPTS_H1 GETPTS_H2 line_points

%[x,y] = getcurpt(GETPTS_AX);
pos = get(GETPTS_FIG, 'CurrentPoint');
x = pos(1);
y = pos(2);
view = getMLRView;
% then use select3d which is slooow, but accurate
hobj = get(get(view.figure, 'CurrentAxes'),'Children');
% make sure we are using the correct object (should be the 3D
% brain). Use end here because with searchForVoxel we plot a
% point on the image and the brain object always seems to be
% last. But if this is not always the case, then we may need
% to do a little more work here to find the correct object
[pos v vi] = select3d(hobj(end));
if isempty(pos),return,end
% convert the index to the coordinates

line_points =[line_points pos];
if(numel(line_points) > 3)
  plot3(line_points(1,:), line_points(2,:), line_points(3,:),'Parent', get(view.figure, 'CurrentAxes'));
end
set([GETPTS_H1 GETPTS_H2], ...
	'XData', x, ...
	'YData', y, ...
	'Visible', 'on');

if (~strcmp(get(GETPTS_FIG, 'SelectionType'), 'normal'))
	% We're done!
	set(GETPTS_H1, 'UserData', 'Completed');
else
	set(GETPTS_FIG, 'WindowButtonDownFcn', 'mlr_getpts_3d(''NextButtonDown'');');
end

%--------------------------------------------------
% Subfunction NextButtonDown
%--------------------------------------------------
function NextButtonDown %#ok

global GETPTS_FIG GETPTS_AX GETPTS_H1 GETPTS_H2 line_points

selectionType = get(GETPTS_FIG, 'SelectionType');
if (~strcmp(selectionType, 'open'))
	% We don't want to add a point on the second click
	% of a double-click
	
	%[newx, newy] = getcurpt(GETPTS_AX);
	pos = get(GETPTS_FIG, 'CurrentPoint');
	newx = pos(1);
	newy = pos(2);
	view = getMLRView;
	% then use select3d which is slooow, but accurate
	hobj = get(get(view.figure, 'CurrentAxes'),'Children');
	% make sure we are using the correct object (should be the 3D
	% brain). Use end here because with searchForVoxel we plot a
	% point on the image and the brain object always seems to be
	% last. But if this is not always the case, then we may need
	% to do a little more work here to find the correct object
	[pos v vi] = select3d(hobj(end));
	% convert the index to the coordinates
	
	if ~isempty(pos)
		line_points =[line_points pos];
		if(numel(line_points) > 3)
			plot3(line_points(1,:), line_points(2,:), line_points(3,:), 'Parent', get(view.figure, 'CurrentAxes'));
		end
	end
	
	
	
	
	x = get(GETPTS_H1, 'XData');
	y = get(GETPTS_H2, 'YData');
	
	set([GETPTS_H1 GETPTS_H2], 'XData', [x newx], ...
		'YData', [y newy]);
	
end

if (~strcmp(get(GETPTS_FIG, 'SelectionType'), 'normal'))
	% We're done!
	set(GETPTS_H1, 'UserData', 'Completed');
end



%----------------------------------------------------
% Subfunction CreatePointer
%----------------------------------------------------
function [pointerShape, pointerHotSpot] = CreatePointer

pointerHotSpot = [8 8];
pointerShape = [ ...
	NaN NaN NaN NaN NaN   1   2 NaN   2   1 NaN NaN NaN NaN NaN NaN
	NaN NaN NaN NaN NaN   1   2 NaN   2   1 NaN NaN NaN NaN NaN NaN
	NaN NaN NaN NaN NaN   1   2 NaN   2   1 NaN NaN NaN NaN NaN NaN
	NaN NaN NaN NaN NaN   1   2 NaN   2   1 NaN NaN NaN NaN NaN NaN
	NaN NaN NaN NaN NaN   1   2 NaN   2   1 NaN NaN NaN NaN NaN NaN
	1   1   1   1   1   1   2 NaN   2   1   1   1   1   1   1   1
	2   2   2   2   2   2   2 NaN   2   2   2   2   2   2   2   2
	NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
	2   2   2   2   2   2   2 NaN   2   2   2   2   2   2   2   2
	1   1   1   1   1   1   2 NaN   2   1   1   1   1   1   1   1
	NaN NaN NaN NaN NaN   1   2 NaN   2   1 NaN NaN NaN NaN NaN NaN
	NaN NaN NaN NaN NaN   1   2 NaN   2   1 NaN NaN NaN NaN NaN NaN
	NaN NaN NaN NaN NaN   1   2 NaN   2   1 NaN NaN NaN NaN NaN NaN
	NaN NaN NaN NaN NaN   1   2 NaN   2   1 NaN NaN NaN NaN NaN NaN
	NaN NaN NaN NaN NaN   1   2 NaN   2   1 NaN NaN NaN NaN NaN NaN
	NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];


