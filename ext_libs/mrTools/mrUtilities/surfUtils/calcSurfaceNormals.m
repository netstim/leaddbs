% calcNormals.m
%
%        $Id$ 
%      usage: normals = calcNormals(surf)
%         by: justin gardner, vectorized for speed by julien besle (2009/09/24)
%       date: 08/10/08
%    purpose: Calculates vertex normals from a surf returned by loadSurfOFF
%             This is computes by averaging the surface normals of each triangle
%             a vertex is a part of
function vertexNormals = calcSurfaceNormals(surf)

% check arguments
if ~any(nargin == [1])
  help calcNormals
  return
end

% first compute the normals to each triangle face.
% this is done with the cross product of two edge vectors

% % % %---------- LOOP VERSION
% % % triNormals = zeros(surf.Ntris,3);
% % % disppercent(-inf,'(calcSurfaceNormal) Computing triangle normals');
% % % for iTri = 1:surf.Ntris
% % %   % get the three vertices of this triangle
% % %   vertex1 = surf.vtcs(surf.tris(iTri,1),:);
% % %   vertex2 = surf.vtcs(surf.tris(iTri,2),:);
% % %   vertex3 = surf.vtcs(surf.tris(iTri,3),:);
% % %   % and compute the surface normal using the cross product
% % %   triNormals(iTri,:) = cross(vertex2-vertex1,vertex2-vertex3);
% % %   triNormals(iTri,:) = triNormals(iTri,:)/norm(triNormals(iTri,:));
% % %   disppercent(iTri/surf.Ntris);
% % % end
% % % disppercent(inf);

%----------- VECTORIZED VERSION (much faster)
% get 2 vector sides for all triangles
vertex2minus1 = surf.vtcs(surf.tris(:,2),:) - surf.vtcs(surf.tris(:,1),:);
vertex2minus3 = surf.vtcs(surf.tris(:,2),:) - surf.vtcs(surf.tris(:,3),:);
% and compute the surface normals using the cross product
triNormals = cross(vertex2minus1,vertex2minus3);
triNormals = triNormals./repmat(sqrt(triNormals(:,1).^2+triNormals(:,2).^2+triNormals(:,3).^2),1,3);
clear vertex2minus1 vertex2minus3


% now compute the normal to each vertex as the average
% normal of all the triangle faces it belongs to

%----------- VECTORIZED VERSION (much faster)
vertexNormals = zeros(surf.Nvtcs,3);   %this will be used to add the normal vectors of all triangles of a given vertex
numberTriangles = zeros(surf.Nvtcs,1); %this will keep track of the number of triangle per vertex
%get all the triangles each vertex belongs to
for iVtx = 1:3; % this loop over the 3 vertices of the triangle matrix
              % could probably be eliminated but requires more bookkeeping of indices
  vertices = surf.tris(:,iVtx);
  while nnz(vertices)
    % find (last) unique occurence of any vertex present in the current column of the triangles matrix
    [whichVertex,whichTriangle] = unique(vertices);
    % remove vertex 0 (which were added to indicate that the vertex has already been found)
    whichTriangle(~whichVertex)=[];
    whichVertex(~whichVertex)=[];
    % add the triangle normals found for these vertices
    vertexNormals(whichVertex,:) = vertexNormals(whichVertex,:) + triNormals(whichTriangle,:);
    numberTriangles(whichVertex) = numberTriangles(whichVertex)+1;
    % in the current column of the triangles matrix, replace the vertices that have just been found by 0
    vertices(whichTriangle)=0;
  end
end
%average triangle normals
vertexNormals = vertexNormals./repmat(numberTriangles,1,3);

    
% % % %---------- LOOP VERSION
% % % disppercent(-inf,'(calcSurfaceNormal) Computing vertex normals');
% % % vertexNormals = zeros(surf.Nvtcs,3);
% % % for iVtx = 1:surf.Nvtcs
% % %   % get which triangles this vertex belongs to
% % %   [triNums edgeNums] = find(iVtx == surf.tris);
% % %   % and then get the mean of the normals to those triangls
% % %   vertexNormals(iVtx,:) = mean(triNormals(triNums,:));
% % %   disppercent(iVtx/surf.Nvtcs);
% % % end
% % % disppercent(inf);

