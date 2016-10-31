% calcDist.m
%
%	$Id$	
%      usage: calcDist(v,method,[coords])
%         by: eli merriam
%       date: 11/28/07
%    purpose: Calculate surface/euclidean distance between points. Updated
%    by Johan to add support for batch processing with 'coords' input,
%    euclidean distance output for 'roi' method.
%
function [dijkstraDistance, euclideanDistance] = calcDist(view,method,coords)
  
% check arguments
if ~any(nargin == [1 2 3])
  help calcDist
  return
end

euclideanDistance = [];
dijkstraDistance = [];
if ieNotDefined('method'); method = 'pairs'; end

% get the base CoordMap for the current flat patch
corticalDepth = mean(viewGet(view, 'corticalDepth'));
baseCoordMap = viewGet(view,'baseCoordMap');
baseHdr = viewGet(view, 'basehdr');

if isempty(baseCoordMap)
  mrWarnDlg('(calcDist) You cannot use this function unless you are viewing a flatpatch with a baseCoordMap');
  return;
end

% use GUI to select points on surface and convert to base coords (default
% behavior)
if ieNotDefined('coords')
  % baseCoords contains the mapping from pixels in the displayed slice to
  % voxels in the current base volume.
  baseCoords = viewGet(view,'cursliceBaseCoords');
  if isempty(baseCoords)
    mrWarnDlg('Load base anatomy before drawing an ROI');
    return;
  end

  disp('(calcDist) Draw line segments for which you want to compute the distance of. Click to add a point. Double-click to end making line segments');
  if strcmp(method, 'roi')
    disp(sprintf('(calcDist) ROI method not implemented'));
  else
    % Select main axes of view figure for user input
    fig = viewGet(view,'figNum');
    gui = guidata(fig);
    set(fig,'CurrentAxes',gui.axis);
    
    % pick some points
    [xi yi] = getpts;

    % draw the lines temporarily
    switch lower(method)
      case {'segments'}
        line(xi,yi, 'linewidth', 3);
      case {'pairs'}
        % ignore the last point if there are an odd number of inputs
        if ~iseven(length(xi))
          xi = xi(1:end-1);
          yi = yi(1:end-1);
        end
        for p=1:2:length(xi)
          line(xi(p:p+1), yi(p:p+1));
        end
    end
    drawnow;  

    % Extract coordinates in base reference frame
    baseX = baseCoords(:,:,1);
    baseY = baseCoords(:,:,2);
    baseZ = baseCoords(:,:,3);
    lineInd = sub2ind(size(baseX), round(yi), round(xi));
    x = baseX(lineInd);
    y = baseY(lineInd);
    z = baseZ(lineInd);
    coords = [x y z];
  end
end

disp('(calcDist) Computing distance');

% load the appropriate surface files
disp(sprintf('Loading %s', baseCoordMap.innerCoordsFileName));
surf.inner = loadSurfOFF(fullfile(baseCoordMap.path, baseCoordMap.innerCoordsFileName));
if isempty(surf.inner)
  mrWarnDlg(sprintf('(calcDist) Could not find surface file %s',fullfile(baseCoordMap.path, baseCoordMap.innerCoordsFileName)));
  return
end
surf.inner = xformSurfaceWorld2Array(surf.inner, baseHdr);

disp(sprintf('Loading %s', baseCoordMap.outerCoordsFileName));
surf.outer = loadSurfOFF(fullfile(baseCoordMap.path, baseCoordMap.outerCoordsFileName));
if isempty(surf.outer)
  mrWarnDlg(sprintf('(calcDist) Could not find surface file %s',fullfile(baseCoordMap.path, baseCoordMap.outerCoordsFileName)));
  return
end
surf.outer = xformSurfaceWorld2Array(surf.outer, baseHdr);

% build up a mrMesh-style structure, taking account of the current corticalDepth
mesh.vertices = surf.inner.vtcs+corticalDepth*(surf.outer.vtcs-surf.inner.vtcs);
mesh.faceIndexList  = surf.inner.tris;

% calculate the connection matrix
[mesh.uniqueVertices,mesh.vertsToUnique,mesh.UniqueToVerts] = unique(mesh.vertices,'rows'); 
mesh.uniqueFaceIndexList = findUniqueFaceIndexList(mesh); 
mesh.connectionMatrix = findConnectionMatrix(mesh);

% get the coordinates of the vertices that are closest to the selected points
[nearestVtcs, distances] = assignToNearest(mesh.uniqueVertices, coords);

% complain if any of the selected points are far away
for i=1:length(distances)
  if (distances>1)
    disp(sprintf('(calcDist) Point %i is %f from the nearest surface vertex', i, distances(i)));
  end
end

% find the distance of all vertices from their neighbours 
scaleFactor = [1 1 1];
D = find3DNeighbourDists(mesh,scaleFactor); 

dijkstraDistance = [];
switch lower(method)
  case {'segments'}
    % calculate length of each segment (1-2, 2-3, 3-4)
    for p=1:length(nearestVtcs)-1
      dist = dijkstra(D, nearestVtcs(p))';
      dijkstraDistance(p) = dist(nearestVtcs(p+1));
      euclideanDistance(p) = norm(diff(coords(p:p+1,:)));
    end 
 case {'pairs'}
    % calculate lengths of a bunch of pairs (1-2, 3-4, 5-6)
    for p=1:2:length(nearestVtcs)
      dist = dijkstra(D, nearestVtcs(p))';
      dijkstraDistance(end+1) = dist(nearestVtcs(p+1));
      euclideanDistance(end+1) = norm(diff(coords(p:p+1,:)));
    end
  case {'roi'}
    % calculate lengths of each point from the first coordinate
    dist = dijkstra(D, nearestVtcs(1))';
    dijkstraDistance = dist(nearestVtcs);
    % JC - euclidean distance between first coord and all others
    normcoords = bsxfun(@minus,coords,coords(1,:));
    euclideanDistance = sqrt(sum(normcoords.^2,2));
  otherwise
    mesh.dist = dijkstra(D, nearestVtcs(1))';
end

% if nargout == 1
%   dijkstraDistance = dist(nearestVtcs);
% end


% $$$ patch('vertices', mesh.uniqueVertices, 'faces', mesh.uniqueFaceIndexList,'FaceVertexCData', mesh.dist,'facecolor', 'interp','edgecolor', 'none');

return;


% patch('vertices', mesh.uniqueVertices, 'faces', mesh.uniqueFaceIndexList, ...
%       'FaceVertexCData', mesh.dist, ...
%       'facecolor', 'interp','edgecolor', 'none');


% baseX = zeros(1,length(xi));
% baseY = zeros(1,length(xi));
% baseZ = zeros(1,length(xi));
% for p=1:length(xi)
%   baseX(p) = baseCoordMap.innerCoords(xi(p), yi(p), 1, 1);
%   baseY(p) = baseCoordMap.innerCoords(xi(p), yi(p), 1, 2);
%   baseZ(p) = baseCoordMap.innerCoords(xi(p), yi(p), 1, 3);
% end

% coords = [baseX' baseY' baseZ'];



