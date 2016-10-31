% mlrGetptsSurface.m
%
%        $Id:$ 
%      usage: [vertices,x,y,z] = mlrGetptsSurface(v,h)
%         by: justin gardner
%       date: 07/11/15
%    purpose: Get user selected points on a surface h
%             works like getpts. Click for points, double-click to end
%             v is the view
%
function [vertices x y z] = mlrGetptsSurface(v,h)

x = [];y = [];z = [];vertices = [];
% check arguments
if ~any(nargin == [2])
  help mlrGetptsSurface
  return
end

% get old user data (used by other interrogator functions)
oldUserData = get(h,'UserData');

% clear the user data
userData = [];
userData.f = viewGet(v,'fignum');
userData.baseSurface = viewGet(v,'baseSurface');
userData.x = [];
userData.y = [];
userData.z = [];
userData.vertices = [];
userData.n = 0;
userData.h = [];
userData.hPath = [];
set(h,'UserData',userData);

% set the cursor
oldPointer = get(userData.f,'Pointer');
[pointerShape, pointerHotSpot] = CreatePointer;
set(userData.f, 'Pointer', 'custom', 'PointerShapeCData', pointerShape, 'PointerShapeHotSpot', pointerHotSpot);

% set the function handler on the surface
oldButtonDown = get(h,'ButtonDownFcn');
set(h,'ButtonDownFcn',@mlrGetptsSurfaceButtonDown);
uiwait(userData.f);

% reset pointer
set(userData.f, 'Pointer', oldPointer);

% get the points from the surface
userData = get(h,'UserData');
x = userData.x;
y = userData.y;
z = userData.z;
vertices = userData.vertices;

% unhook the button down handler
set(h,'UserData',oldUserData);
set(h,'ButtonDownFcn',oldButtonDown);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   mlrGetptsSurfaceButtonDown   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrGetptsSurfaceButtonDown(hObject,e,handles,varargin)

% get the points already clicked
userData = get(hObject,'UserData');
if ~isfield(userData,'hPath'),return,end

% remember clicked point
userData.lastClick = e.IntersectionPoint;
% get the vertex it corresponds to
% compute distance to every vertex and pick the vertex that is closest to the intersectoin point
[minDist userData.vertices(end+1)] = min(sum((e.Source.Vertices-repmat(e.IntersectionPoint',1,size(e.Source.Vertices,1))').^2,2));
% save the x,y and z of this point
userData.x(end+1) = e.Source.Vertices(userData.vertices(end),1);
userData.y(end+1) = e.Source.Vertices(userData.vertices(end),2);
userData.z(end+1) = e.Source.Vertices(userData.vertices(end),3);
% update count
userData.n = userData.n+1;
% remove last drawing of points
delete(userData.h);
delete(userData.hPath);
% save data in object
set(hObject,'UserData',userData);

% check if we hit the same point twice (cheap double-click check - since
% the surface object does not seem to have a proper double-click Selection-Type field
if (userData.n && (userData.x(end) == e.IntersectionPoint(1)) && (userData.y(end) == e.IntersectionPoint(2))) || isequal(get(userData.f,'SelectionType'),'alt')
  % double-click, we got all points
  if userData.n > 1
    % plot the points
    userData.h = plot3(userData.x([1:end-1 1]),userData.y([1:end-1 1]),userData.z([1:end-1 1]),'wo','Parent', get(userData.f, 'CurrentAxes'),'MarkerSize',6,'Color','w','MarkerFaceColor','w');
    % plot the path between this point and the first point
    pathVertices = mlrGetPathBetween(userData.baseSurface,userData.vertices(end-1),userData.vertices(1));
    userData.hPath = plot3(userData.baseSurface.vtcs(pathVertices,1),userData.baseSurface.vtcs(pathVertices,2),userData.baseSurface.vtcs(pathVertices,3),'w-');
    drawnow;
  end
  % resume above function
  uiresume(userData.f);
else
  % plot the points
  userData.h = plot3(userData.x,userData.y,userData.z,'wo','Parent', get(userData.f, 'CurrentAxes'),'MarkerSize',6,'Color','w','MarkerFaceColor','w');
  if userData.n > 1
    % plot the path between this point and the last point
    pathVertices = mlrGetPathBetween(userData.baseSurface,userData.vertices(end-1),userData.vertices(end));
    userData.hPath = plot3(userData.baseSurface.vtcs(pathVertices,1),userData.baseSurface.vtcs(pathVertices,2),userData.baseSurface.vtcs(pathVertices,3),'w-');
    drawnow;
  end
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
