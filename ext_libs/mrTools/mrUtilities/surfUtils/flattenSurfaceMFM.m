function [unfoldMeshSummary] = flattenSurfaceMFM(mesh, startCoords, perimDist);
% Author: Wade
% Purpose:
%    Unfolds a properly triangulated mesh using Floater/Tutte's sparse
%    linear matrix method. Then maps layer one grey nodes to the flattened
%    mesh and squeezes remaining grey layers to L1 to generate a flat.mat
%    file for use with mrLoadRet/mrVista
%
%    1. Read in the white matter boundary mesh created by mrGray.  The data
%    in the white matter mesh are organized into the mesh data structure.
%    2. The second section pulls out the portion of this entire mesh that will
%    be unfolded, creating the unfoldMesh structure.  This is the portion of 
%    the full mesh that is within the criterion distance from the startCoords.  
%    3. The data in the unfoldMesh are flattened.
%    4. The flat positions are adjusted to make the spacing more nearly
%    like the true distances.
%    5. Gray matter data to the flat map positions, building gLocs2d and
%    gLocs3d that are used by mrVista/mrLoadRet.
%
% See Also:  mrFlatMesh (the GUI that calls this).
%
% Stanford University

NPERIMS            = 1;
SAVE_INTERMEDIATE  = 0;
NUMBLURSTEPS       = 1; % Blur the curvature map a little in a connection-dependent manner
scaleFactor = [1 1 1];
statusHandle = 0;
busyHandle = 0;
showFigures = 0;

if (length(startCoords)~=3), error ('Error: you must enter 3 start coords'); end
if ieNotDefined('spacingMethod'), spacingMethod = 'None'; end
if ieNotDefined('gridSpacing'), gridSpacing = 0.4; end
%verbose = strcmp(mrGetPref('Verbose'),'Yes');
verbose = 1;

[adjustSpacing, spacingMethod] = xlateSpacing(spacingMethod);

startTime=now;
str = sprintf('****** mrFlatMesh %s *****',datestr(now));
if verbose,disp(sprintf(str));end

scaleFactorFromMesh = [1 1 1];

if ieNotDefined('scaleFactorFromMesh')
  errordlg('No mesh scale factor. You must create a new (modern) mesh.'); 
end
str = sprintf('Scale factor %.03f,%.03f,%.03f',scaleFactorFromMesh);
if verbose,disp(sprintf(str));end
scaleFactor=scaleFactorFromMesh; % mm per voxel

%-------------------------------------------------
% Step 1.  We build the mesh from the mrGray output mesh.  This contains
% all of the segmentation information.  This should be a separate routine.
%-------------------------------------------------

% The entries in mesh.vertices aren't unique. That is, the same
% vertex can be represented by two different entries.  This happens
% because strips always contain shared vertices.
%
% This next series of functions finds the unique vertices and
% and builds hash tables for vert->uniqueVerts and uniqueVerts->vert.
% Using these tables, we can adjust the list of colors and faces
% to represent each vertex by a unique index.

% Get rid of identical points
[mesh.uniqueVertices,mesh.vertsToUnique,mesh.UniqueToVerts]= unique(mesh.vertices,'rows'); 

% try to assign the colors that came in from mrGray
mesh.uniqueCols=mesh.rgba(mesh.vertsToUnique,:); 

% this gets rid of faces with multiple duplicate vertices and different permutations of the same face
mesh.uniqueFaceIndexList=findUniqueFaceIndexList(mesh); 

% Now we find the connection matrix: 
% a sparse matrix of nxn points where M(i,j) is 1 if i and j are connected
if verbose,statusStringAdd(statusHandle,'Finding connection matrix.');end
[mesh.connectionMatrix]=findConnectionMatrix(mesh);
sumCon=sum(mesh.connectionMatrix);
mesh.uniqueCols=mesh.uniqueCols(:,1);

if(NUMBLURSTEPS>0)
  if verbose,statusStringAdd(statusHandle,'Blurring colormap.');end
  for t=1:NUMBLURSTEPS
    mesh.uniqueCols=(mesh.connectionMatrix*mesh.uniqueCols(:))./sumCon(:);
  end
end

str=sprintf('%d connections found',length(mesh.connectionMatrix));
if verbose,disp(sprintf(str));end

% At this point, we can use the connection matrix to perform some fast
% smoothing on the curvature map. Also possibly a relaxation / smoothing on
% the actual mesh?
for t=1:NUMBLURSTEPS
  mesh.uniqueCols = connectionBasedSmooth(mesh.connectionMatrix, mesh.uniqueCols(:,1));
end

if verbose,statusStringAdd(statusHandle,'Checking group perimeter.');end

% Check to make sure that this is a clean surface: no edge points yet.
edgeList=findGroupPerimeter(mesh,1:length(mesh.uniqueVertices));
if (~isempty(edgeList))
  % it turns out that some surfaces can not be open and still not cause a problem
  % especially if the edge points are not within the location being flattened. So,
  % making this a warning, rather than a full on error (perhaps a better test
  % would be to check whether the vertices in question live witin the are to be
  % flattened, but since this seems to be a rather rare occurrence (happened for
  % us with a freesurfer surface with a bright spot in the orbital frontal cortex).
  %error('Error - initial mesh is not closed!');
  %return;
  dispHeader('(flattenSurfaceMFM) ERROR - initial mesh is not closed!',80);
  disp(sprintf('(flattenSurfaceMFM) If this happens within the region to be flattened, this may cause a problem'));
  disp(sprintf('(flattenSurfaceMFM) But otherwise, your flat surface may be ok. The following is a list of the'));
  disp(sprintf('(flattenSurfaceMFM) coordinates in the 3D canonical volume from which the problem seems to be occurring'));
  disp(sprintf('(flattenSurfaceMFM) You are ENCOURAGED to go check these points (load the canonical in mlrVol and look'));
  disp(sprintf('(flattenSurfaceMFM) at the matching X,Y,Z locations and make sure this is out of the region that you want'));
  disp(sprintf('(flattenSurfaceMFM) to flatten. You may also consider trying to fix the problem (perhaps scan another canonical'));
  disp(sprintf('(flattenSurfaceMFM) and or change parameters of the segmentation software to get surfaces without this problem'));
  for i = 1:length(edgeList(:))
    disp(sprintf('(flattenSurfaceMFM) (X,Y,Z): (%i, %i, %i)',round(mesh.uniqueVertices(edgeList(i),1)),round(mesh.uniqueVertices(edgeList(i),2)),round(mesh.uniqueVertices(edgeList(i),3))));
  end
  dispHeader('',80);
else
  str = sprintf('Initial mesh is closed.');
  if verbose,disp(sprintf(str));end
end

if verbose,statusStringAdd(statusHandle,'Finding closest mesh point.');end

% Find the nearest mesh point to the startCoords (Euclidian distance).
[startNode,snDist]=assignToNearest(mesh.uniqueVertices,startCoords); 

% Print the distance from the gray matter and warn if you're more than 15
% voxels away
str=sprintf('Start node %d selected at %d voxel units from input coords.',startNode,sqrt(snDist));
if verbose,disp(sprintf(str));end
if (sqrt(snDist)>15)
  beep;
  str=sprintf('** Warning: mesh node far from start coord. Expect trouble.');
  if verbose,statusStringAdd(statusHandle,str);end
end

% Find the distance of all nodes from the start node so that we can unfold
% just a sub-region of the whole mesh 
if verbose,statusStringAdd(statusHandle,'Finding distances from start node');end

% To this point, we are in a voxel framework: 
% Everything has been scaled in the mrReadMrM function. 
%
% D is the connection matrix using the true (non-squared) distance.
D = sqrt(find3DNeighbourDists(mesh,scaleFactor)); 

str=sprintf('Mean inter-node distance: %.03f',full(sum(sum(D)))/nnz(D));
if verbose,disp(sprintf(str));end

% Find distances from the startNode to all the nodes
% Have replaced mrManDist with 'dijkstra' mex file to get around potential rounding errors.
mesh.dist=dijkstra(D,startNode);

% We now have the basic mesh information.  We are starting to identify the
% perimeter and inside nodes for flattening. We generate a perimeter, based
% on the user's choice of perimeter distasnce  by thresholding these
% distances
mesh.perimDist=perimDist;
messageString=sprintf('Perimeter minimum distance=%d',perimDist);

if verbose,statusStringAdd(statusHandle,messageString);end
if verbose,statusStringAdd(statusHandle,'Using threshold to find perimeter(s).');end
if verbose,statusStringAdd(statusHandle,'Defining perimeters.');end

[perimeterEdges,eulerCondition]=findLegalPerimeters(mesh,perimDist);

% The routine above can generate islands - fix it by zeroing the connection matrix for the largest perimeter
% and then doing a flood fill with no limits to generate the inside group

% FIND ALL THE PERIMETERS AND TAKE THE BIGGEST ONE
messageString=sprintf('Euler number for this set=%d',eulerCondition);
if verbose,statusStringAdd(statusHandle,messageString);end

uniquePerimPoints=unique(perimeterEdges(:));
messageString=sprintf('%d unique perimeter points found.',length(uniquePerimPoints));
if verbose,statusStringAdd(statusHandle,messageString);end

[orderedUniquePerimeterPoints,biggest]=orderMeshPerimeterPointsAll(mesh,perimeterEdges);
nPerims=size(orderedUniquePerimeterPoints);
orderedUniquePerimeterPoints=orderedUniquePerimeterPoints{biggest}.points;

% DO THE FLOOD FILL TO FILL UP THE INNER HOLES
tempConMat=mesh.connectionMatrix; % save it for later
mesh.connectionMatrix(orderedUniquePerimeterPoints,:)=0;
mesh.connectionMatrix(:,orderedUniquePerimeterPoints)=0;
if verbose,statusStringAdd(statusHandle,'Doing flood fill.');end

[insideNodes,insideNodeStruct]=floodFillFindPerim(mesh,Inf,startNode,busyHandle);
insideNodes=[insideNodes(:);orderedUniquePerimeterPoints(:)];
mesh.connectionMatrix=tempConMat;
[perimeterEdges,eulerCondition]=findGroupPerimeter(mesh,insideNodes);
str = sprintf('Euler condition=%d',eulerCondition);
if verbose,disp(sprintf(str));end

% We now have a fully-connected mesh, and we have identified perimeter and
% inside nodes.  We are ready to build the portion of the mesh that we will
% unfold.
%

%-------------------------------------------------
% Step 2.  We extract the part of the mesh that we will unfold.  This part
% is defined by the distance and start node selected by the user.
% We now have a fully-connected mesh, and we have identified perimeter and
% inside nodes.  We are ready to build the portion of the mesh that we will
% unfold.
%-------------------------------------------------
[unfoldMesh, nFaces] = mfmBuildSubMesh(mesh, perimeterEdges, insideNodes, ...
                                       orderedUniquePerimeterPoints, statusHandle,busyHandle);


%-------------------------------------------------
% Step 3.  Unfold the unfoldMesh.  This is the key mathematical step in the
% process.
%-------------------------------------------------

verbose = strcmp(mrGetPref('Verbose'),'Yes');

% Find the N and P connection matrices
if verbose,statusStringAdd(statusHandle,'Finding sub-mesh connection matrix.');end
[N, P, unfoldMesh.internalNodes] = findNPConnection(unfoldMesh);

% Here we find the squared 3D distance from each point to its neighbours.
[unfoldMesh.distSQ]=find3DNeighbourDists(unfoldMesh,scaleFactor);   

fullConMatScaled=scaleConnectionMatrixToDist(sqrt(unfoldMesh.distSQ));

% Now split the full conMat up until N and P
N = fullConMatScaled(unfoldMesh.internalNodes,unfoldMesh.internalNodes);
P = fullConMatScaled(unfoldMesh.internalNodes,unfoldMesh.orderedUniquePerimeterPoints);

% Assign the initial perimeter points - they're going to go in a circle for now...
numPerimPoints=length(unfoldMesh.orderedUniquePerimeterPoints);

if verbose,statusStringAdd(statusHandle,'Assigning perimeter points');end

% Can set distances around circle to match actual distances from the center. 
unfoldMesh.X_zero = assignPerimeterPositions(perimDist,unfoldMesh); 


% THIS IS WHERE WE SOLVE THE POSITION EQUATION - THIS EQUATION IS THE HEART OF THE ROUTINE!
if verbose,statusStringAdd(statusHandle, 'Solving position equation.');end
X =(speye(size(N)) - N) \ (sparse(P * unfoldMesh.X_zero));

% Remember what these variables are: 
% X: 2d locations of internal points
% X_zero : 2d locations of perimeter
% N sparse connection matrix for internal points
% P sparse connection matrix between perimeter and internal points
% (Note - the connection matrix between the perimeter points is implicit in
% their order - they are connected in a ring)
% The 3D coordinates o

unfoldMesh.N=N;
unfoldMesh.P=P;
unfoldMesh.X=X;

% Find the face areas for the unfolded and folded versions of the unfoldMesh
% This is a good error metric. We'll save this out with the flat.mat file.
% (ARW)
%
% Do this by calling findFaceArea. We call it twice, once with the 3D vertices and once with 
% a pseudo-3D vertex set with the 3rd dimension set to 0
% The areaList3D and errorList are saved out,
% but I don't know where they are used.  Perhaps Bob? (BW)

statusStringAdd(statusHandle,['Calculating face area distortions']);

% This seems like an important piece of code, the ordering used to define
% unfolded2D is complicated.  We need comments and it would be better to
% have it in a function.
unfolded3D = unfoldMesh.uniqueVertices;
indices2D  = unfoldMesh.internalNodes;
unfolded2D(indices2D,1:2) = full(unfoldMesh.X);
indices2D  = unfoldMesh.orderedUniquePerimeterPoints;
unfolded2D(indices2D,1:2) = full(unfoldMesh.X_zero);
unfolded2D(:,3) = 0;

% Is this related to RFD stuff?
areaList3D=findFaceArea(unfoldMesh.connectionMatrix,unfolded3D,unfoldMesh.uniqueFaceIndexList);
areaList2D=findFaceArea(unfoldMesh.connectionMatrix,unfolded2D,unfoldMesh.uniqueFaceIndexList);

% Hmm.  We get a divide by zero somtimes, indicating that the 2D area is
% zero.  That can't be good.  I protected this by adding eps.  But that is
% just to spare the user.
errorList = areaList3D./(areaList2D + eps);
zeroAreaList = find(areaList2D == 0);
if ~isempty(zeroAreaList), fprintf('Zero 2D area (nodes): %.0f',zeroAreaList); end

% In order to plot this as a nice picture, we want to find the (x,y) center
% of mass of each 2D face.  
% I think cogs means center-of-gravity (BW).  I am not sure we use this
% any more, and I am not sure we use the areaList stuff above, either.
if (showFigures)
  figure; subplot(1,3,1)
  [areaErrorMap,meanX,meanY] =  mfmAreaErrorMap(unfoldMesh, nFaces, unfolded2D,errorList);
end

%--------------------------------------------
%Step 4.  The unfoldMesh is complete. Now we adjust the spacing of the
%points so they don't bunch up too much.  The method is to find a cartesian
%grid within the data, find a Delaunay triangulation of this grid, and then
%use a series of affine transformations to transform each of the triangles
%to an equal area representation with the proper grid topology.  This is
%explained in more detail in flatAdjustSpacing
%--------------------------------------------

if (adjustSpacing)
  str = sprintf('Spacing points (method = %s)',spacingMethod);
  if verbose,disp(sprintf(str));end
  
  unfoldMesh.maxFractionDist = 0.8;
  unfoldMesh.gridSpacing = gridSpacing;
  unfoldMesh.locs2d = unfolded2D(:,1:2);
  unfoldMesh.startCoords = startCoords;
  unfoldMesh.scaleFactor = scaleFactor;
  
  % The user can select Cartesian or Polar spacing methods.  Only
  % Cartesian is working now, though, I think -- BW
  [newLocs2d,goodIdx] = flatAdjustSpacing(unfoldMesh,spacingMethod);
end

if (showFigures)
  statusStringAdd(statusHandle,['Displaying unfolded mesh']);
  unfoldMeshFigure(unfoldMesh); 
end

% Finally - the mapping of grey to mesh points takes place using the entire mesh. 
% Therefore, we need to generate X for the mesh as well as the unfold mesh

% insideNodes is an array of indices into the original (non-bounded) mesh.
% Each entry in insideNodes relates a point in the unfoldMesh to a point in
% the original mesh 

% Recover the perimeter and internal points
mesh.X=zeros(length(mesh.uniqueVertices),2);

% In order to deal with a cropped mesh (as generated by flatAdjustSpacing)
% we need to compute the ... Alex?
if (adjustSpacing)
  
  mesh.X(insideNodes(goodIdx),:)=newLocs2d;
  hasCoords=[insideNodes(goodIdx)];
  
else
  
  unfoldToOrigPerimeter=insideNodes(unfoldMesh.orderedUniquePerimeterPoints);
  unfoldToOrigInside=insideNodes(unfoldMesh.internalNodes);
  
  mesh.X(unfoldToOrigPerimeter,:)=unfoldMesh.X_zero;
  mesh.X(unfoldToOrigInside,:)=unfoldMesh.X;
  hasCoords=[unfoldToOrigPerimeter(:);unfoldToOrigInside(:)];
  
end

coords=mesh.X(hasCoords,:);
dists=mesh.dist(hasCoords);

% use griddata to image the distance map
warning off MATLAB:griddata:DuplicateDataPoints;
mesh.distMap=makeMeshImage(coords,dists,128);
warning on MATLAB:griddata:DuplicateDataPoints;

ZI=mesh.distMap;

if (showFigures), unfoldDistFigure(mesh); end

% Record which nodes in the big mesh are in the unfold
mesh.insideNodes=insideNodes;


endTime=now;

%----------------------------------
% Step 6.
% Save stuff
%----------------------------------

unfoldMeshSummary.startCoords=startCoords;
unfoldMeshSummary.connectionMatrix = unfoldMesh.connectionMatrix;
unfoldMeshSummary.uniqueVertices = unfoldMesh.uniqueVertices;
unfoldMeshSummary.uniqueFaceIndexList = unfoldMesh.uniqueFaceIndexList;
unfoldMeshSummary.internalNodes = unfoldMesh.internalNodes;
unfoldMeshSummary.orderedUniquePerimeterPoints = unfoldMesh.orderedUniquePerimeterPoints;
unfoldMeshSummary.scaleFactor = scaleFactor;
unfoldMeshSummary.locs2d = unfolded2D(:,1:2);

unfoldMeshSummary.insideNodes = mesh.insideNodes;
unfoldMeshSummary.UniqueToVerts = mesh.UniqueToVerts;
unfoldMeshSummary.vertsToUnique = mesh.vertsToUnique;

if showFigures
  unfoldPlotFlatMesh(unfoldMeshSummary.locs2d, unfoldMeshSummary.uniqueFaceIndexList, unfoldMesh.uniqueCols)
end

str = sprintf('****** End mrFlatMesh  %s *****',datestr(endTime));
if verbose,disp(sprintf(str));end

return;

%----------------------------------
function [adjustSpacing, spacingMethod] = xlateSpacing(spacingMethod)
switch(spacingMethod)
  case 'None'
    adjustSpacing = 0;
  case 'Cartesian'
    adjustSpacing = 1;
    spacingMethod = 'cartesian';
  case 'Polar'
    warndlg('Polar not yet implemented.  Using Cartesian equal.');
    adjustSpacing = 1;
    spacingMethod = 'polar';
end
return;

%----------------------------------
function unfoldMeshFigure(unfoldMesh);
subplot(1,3,1)
hold off;
gplot(unfoldMesh.N,unfoldMesh.X);
title ('Unfolded mesh'); axis equal; axis off; zoom on;
axis image
return;

%----------------------------------
function unfoldDistFigure(mesh)
subplot(1,3,2)
imagesc(mesh.distMap);
axis image; colormap hot; title('Manifold distance map'); colorbar('South');

% Added the abs to cope with the fact that we sometimes get -ve numbers here...
% imagesc(abs(log(areaErrorMap))); 
% axis image;
% colorbar;
% title('Log 3D/2D area');
% colormap hot;

return;


% %----------------------------------
% function unfoldPlotCurvatureMap(gLocs2d,numNodes,mRange,layerCurvature)
% subplot(1,3,2)
% [y x]=meshgrid(mRange);

% warning off MATLAB:griddata:DuplicateDataPoints;
% fl = griddata(gLocs2d(1:numNodes(1),1),gLocs2d(1:numNodes(1),2),layerCurvature{1},x,y);
% warning on MATLAB:griddata:DuplicateDataPoints;

% figure;  
% imagesc((rot90(fl))); colormap gray; title('Curvature map'); axis image

% return;

%----------------------------------
function unfoldPlotFlatMesh(vertices, faces, curv)
subplot(1,3,3)
curv(find(curv>0))  = .15;
curv(find(curv<=0)) = .3;
p(1) = patch('vertices', vertices, 'faces', faces, ...
             'facecolor', 'interp', ...
             'edgecolor', 'none', ...
             'FaceVertexCData', [curv curv curv]);
axis image
axis off
return;
