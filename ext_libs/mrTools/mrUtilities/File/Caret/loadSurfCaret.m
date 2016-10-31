% loadSurfCaret.m
%
%        $Id: loadSurfCaret.m 1281 2008-08-22 00:33:11Z justin $ 
%      usage: surf = loadSurfCaret(coordFilename,<topoFilename>,<'dispSurface=0'>)
%         by: justin gardner
%       date: 08/10/08
%    purpose: load a caret surface, if you want to display, set dispSurface=1
%             returned surface is in same format as returned by loadSurfOFF
%             There are a few options:
%             displaySurface=0 (calculates curvature and displays surface in a figure
%             zeroBased=0 (returns vertex coordinates as 0 based coordinates)
%             coordShift=[x y z] (shift the vertex coordinates by these amounts. Typically
%                for freeSurfer coordinates are shifted by half the volume size)
%             xform=[] (A 4x4 xform matrix that will be applied to the surface coordinates)
%
function surf = loadSurfCaret(coordFilename,topoFilename,varargin)

% check arguments
if ~any(nargin == [1 2 3 4 5 6 7])
  help loadSurfCaret
  return
end

dispSurface = [];
xform = [];
coordShift = [];
zeroBased = [];
getArgs(varargin,{'dispSurface=0','coordShift=[]','xform=[]','zeroBased=0'});

% set extensions
%coordFilename = setext(coordFilename,'coord');

% open coord file
coord = openCaretFile(coordFilename);
if isempty(coord),return,end

% one argument, get topo from coord header
if ieNotDefined('topoFilename')
  if ~isfield(coord,'topo_file')
    topoFilename = setext(coordFilename,'topo');
  else
    topoFilename = coord.topo_file;
  end
end

% open topo file
topo = openCaretFile(topoFilename);
if isempty(topo),return,end
if ~isfield(topo,'num_tiles')
  disp(sprintf('(loadSurfCaret) %s is not a topo file',topoFilename));
  return
end

% create surface
surf.Nvtcs = coord.num_nodes;
surf.Ntris = topo.num_tiles;
surf.Nedges = surf.Nvtcs+surf.Ntris-2;
surf.vtcs = coord.data;
surf.tris = topo.data;

% transform coordinates
if ~isempty(xform)
  vtcs = surf.vtcs';
  vtcs(4,:) = 1;
  vtcs = xform*vtcs;
  surf.vtcs = vtcs(1:3,:)';
end

% zero based coordinates
if zeroBased
  surf.vtcs = surf.vtcs-1;
end

if ~isempty(coordShift)
  if length(coordShift) ~= 3
    disp(sprintf('(loadSurfCaret) Coord shift should be [x y z]'));
    return
  end
  % otherwise shift
  for i = 1:3
    surf.vtcs(:,i) = surf.vtcs(:,i)+coordShift(i);
  end
end

% display the surface
if dispSurface
  surf.m = calcCurvature(surf);
  mlrSmartfig('loadSurfCaret','reuse');clf
  patch('vertices', surf.vtcs, 'faces', surf.tris,'FaceVertexCData', surf.m,'facecolor','interp','edgecolor','none');
  axis equal
end


