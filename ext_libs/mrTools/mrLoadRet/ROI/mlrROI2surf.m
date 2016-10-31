% mlrROI2surf.m
%
%        $Id:$ 
%      usage: surf = mlrROI2surf(roi,roiSurf,<corticalDepth=0.5>)
%         by: justin gardner
%       date: 07/30/15
%    purpose: function that will convert an mlr ROI into a surf (triangulated mesh with
%             vertices and triangles). The roi is the roi you want to convert
%             and roiSurf is the surface you want to base the conversion
%             on (i.e. it will pull out the piece of roiSurf that
%             corresponds to roi)
%
%             v = newView;
%             v = loadROI(v,'l_mt.mat');
%             roi = viewGet(v,'roi','l_mt');
%             v = loadAnat(v,'jg_left_GM');
%             roiSurf = view(v,'base');
%             surf = mlrROI2surf(roi,roiSurf);
%
%             Will default to getting the surface at cortical depth 0f 0.5
%             Can get at different depths by doing
%
%             surf = mlrROI2surf(roi,roiSurf,'corticalDepth=0');
%
function surf = mlrROI2surf(roi,roiSurf,varargin)

% check arguments
if nargin < 2
  help mlrROI2surf
  return
end

% get arguments
getArgs(varargin,{'corticalDepth=0.5'});

if ~isbase(roiSurf) || (roiSurf.type ~= 2)
  disp(sprintf('(mlrROI2surf) Invalid surface passed in'));
  return 
end

% get vertices and some information
vtcs = getSurfVerticesAtCorticalDepth(roiSurf,corticalDepth);
nVertices = size(vtcs,1);
baseDims = roiSurf.coordMap.dims;

% Note that we use rounded coordinates here - this assumes that rois are only
% kept to at most 1x1x1 mm resolution and should change later
vtcs = round(vtcs);
roiCoords = round(roi.coords);

% linearize for faster comparisons
linearVtcs = mrSub2ind(baseDims,vtcs(:,1),vtcs(:,2),vtcs(:,3));
linearRoiCoords = mrSub2ind(baseDims,roiCoords(1,:),roiCoords(2,:),roiCoords(3,:))';

% now check to see which vertices in the base are matched to the roi
% vertices. Note that since we have rounded there can be multiple
% surface vertices that match an roi vertex
vertexIndices = find(ismember(linearVtcs,linearRoiCoords));

matchingTris = roiSurf.coordMap.tris;
% get all triangles associated with these vertices
for i = 1:length(vertexIndices)
  % set any matching vertices to nan
  matchingTris(matchingTris==vertexIndices(i)) = nan;
end
[triNum edgeNum] = find(isnan(matchingTris));
surf.tris = roiSurf.coordMap.tris(triNum,:);

% map into new vertex numbers
unmatchedTris = surf.tris;
for i = 1:length(vertexIndices)
  matchVertices = (surf.tris==vertexIndices(i));
  unmatchedTris(matchVertices) = nan;
end
% add unmatched vertices in
vertexIndices = unique(union(vertexIndices,unmatchedTris(~isnan(unmatchedTris))));
for i = 1:length(vertexIndices)
  matchVertices = (surf.tris==vertexIndices(i));
  surf.tris(matchVertices) = i;
end
% get the vertices
surf.vtcs = vtcs(vertexIndices,:);

% debug code to draw surface and roi surface to visualize what is going on
if 0
  mlrSmartfig('mlrROI2surf','reuse');clf;
  patch('vertices', roiSurf.coordMap.innerVtcs + (roiSurf.coordMap.outerVtcs-roiSurf.coordMap.innerVtcs)*0, 'faces', roiSurf.coordMap.tris,'facecolor','interp','edgecolor','none','FaceAlpha',1,'FaceVertexCData',roiSurf.data');
  hold on
  patch('vertices',surf.vtcs,'faces',surf.tris,'facecolor','interp','FaceVertexCData',repmat([1 0 0],size(surf.vtcs,1),1),'edgecolor','none');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getSurfVerticesAtCorticalDepth    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vtcs = getSurfVerticesAtCorticalDepth(s,corticalDepth)

disp(sprintf('corticalDepth: %f',corticalDepth));
% interpoate coordinates. 
vtcs = s.coordMap.innerVtcs + (s.coordMap.outerVtcs - s.coordMap.innerVtcs) * corticalDepth;
vtcs = squeeze(vtcs);
