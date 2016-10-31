% mlrGetPathBetween.m
%
%        $Id:$ 
%      usage: pathList = mlrGetPathBetween(surf,v1,v2)
%         by: justin gardner
%       date: 07/12/15
%    purpose: Gets shortest path between two vertices given surface
%
%       e.g.: baseSurface = viewGet(v,'baseSurface');
%             mlrGetPathBetween(baseSurface, 100,150);
%
function pathList = mlrGetPathBetween(surf,v1,v2)

% check arguments
if ~any(nargin == [3])
  help mlrGetPathBetween
  return
end

% set up the mesh so that we can get the distanceMatrix
mesh.uniqueVertices = surf.vtcs;
mesh.uniqueFaceIndexList = surf.tris;
mesh.connectionMatrix = findConnectionMatrix(mesh);
distanceMatrix = find3DNeighbourDists(mesh);

% now use dijkstrap to get the shortest path to v1
[d predList] = dijkstrap(distanceMatrix,v1);

% now go backwords from v2 and find the path to v1
pathList = v2;
while pathList(end) ~= v1
  pathList(end+1) = predList(pathList(end));
end
