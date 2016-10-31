% makeFlat.m
%
%       $Id$	
%      usage: makeFlat()
%         by: eli merriam
%       date: 09/27/07
%    purpose: 
%
function retval = makeFlat(view, overlayNum, scan, x, y, s, roi) 


% check arguments
if ~any(nargin == [4 5 6 7])
  help makeFlat
  return
end

% get base info
baseCoordMap = viewGet(view,'baseCoordMap');
baseCoordMapPath = viewGet(view,'baseCoordMapPath');
startPoint = viewGet(view,'mouseDownBaseCoords');
baseType = viewGet(view,'baseType');

% some other variables
defaultRadius = 75;
viewNum = viewGet(view, 'viewNum');

% parse the parameters
paramsInfo = {};
if ~isempty(baseCoordMap)
  % get the inner and outer coordinate surfaces.
  outerCoordsFileName = setext(baseCoordMap.outerCoordsFileName,'off');
  innerCoordsFileName = setext(baseCoordMap.innerCoordsFileName,'off');

  % now see if user has to go find file;
  if ~isfile(fullfile(baseCoordMapPath,innerCoordsFileName))
    startPathStr = baseCoordMapPath;
    filterspec = {'*.off','SurfRelax off file';'*WM*.off', 'SurfRelax OFF WM file'; '*.*','All files'};
    title = 'Choose inner (WM) OFF surface file';    
    innerCoordsFileName = mlrGetPathStrDialog(startPathStr,title,filterspec,'off');
    % Aborted
    if isempty(innerCoordsFileName),return,end
    % get path and name
    [baseCoordMapPath innerCoordsFileName] = fileparts(innerCoordsFileName);
  end
  
  % Make flat structure
  flat.path = baseCoordMapPath;
  flat.parentSurfaceName = innerCoordsFileName;
  flat.startPoint = startPoint;
  flat.radius = defaultRadius;

  % and remember other field settings
  curvFileName = baseCoordMap.curvFileName;
  anatFileName = baseCoordMap.anatFileName;

else
  % Ask the user to choose the inner surface
  startPathStr = mrGetPref('volumeDirectory');
  filterspec = {'*.off','SurfRelax OFF file';'*WM*.off','SurfRelax off gray matter file'; '*.*','All files'};
  title = 'Choose inner (WM) OFF file';
  innerCoordsFileName = mlrGetPathStrDialog(startPathStr,title,filterspec,'off');
  % Aborted
  if isempty(innerCoordsFileName),return,end
  % get path and name
  [flat.path flat.parentSurfaceName] = fileparts(innerCoordsFileName);
  % Make flat structure
  flat.startPoint = startPoint;
  flat.radius = defaultRadius;
  % and set other field settings
  innerCoordsFileName = getLastDir(innerCoordsFileName);
  outerCoordsFileName = [];
  curvFileName = [];
  anatFileName = [];
end

% now bring up mrFlatViewer to set initial round of parameters
params = mrFlatViewer(flat,outerCoordsFileName,innerCoordsFileName,curvFileName,anatFileName,viewNum);
if isempty(params)
  disp('(makeFlat) User pressed Cancel.');
  return;
end

% get some other parameters
paramsInfo = {};
paramsInfo{end+1} = {'flatRes', 2, 'resolution of flat patch', 'round=1', 'minmax=[1 10]', 'incdec=[-1 1]', 'The resolution of the flat patch -- a value of 2 doubles the resolution'};
paramsInfo{end+1} = {'threshold', 1, 'type=checkbox', 'Thresholding the surface makes the background two-tone (binary curvature)'};
paramsInfo{end+1} = {'flattenMethod', {'mrFlatMesh','surfRelax'},'type=popupmenu','Use either surfRelax or mrFlatMesh'};
extraParams = mrParamsDialog(paramsInfo, 'makeFlat', []);
if isempty(extraParams)
  disp('(makeFlat) User pressed Cancel.');
  return;
end

% add extraParams to params
extraParamsFieldNames = fieldnames(extraParams);
for i = 1:length(extraParamsFieldNames)
  params.(extraParamsFieldNames{i}) = extraParams.(extraParamsFieldNames{i});
end

% check for directory
if ~isdir(params.path)
  params.path = uigetdir(mrGetPref('volumeDirectory'),'Find anatomy directory');
  % user hits cancel
  if params.path == 0
    return
  end
end

% now go flatten
if strcmp(extraParams.flattenMethod, 'surfRelax')
  disp(sprintf('Flattening using SurfRelax'));
  flatBase = myCutAndFlatten(params);
elseif strcmp(extraParams.flattenMethod, 'mrFlatMesh')
  disp(sprintf('Flattening using the mrVista mrFlatMesh'));
  flatBase = runMrFlatMesh(params);
end

% make it into a MLR4 base anatomy
if ~isempty(flatBase)
  % make sure name does not have .off on it
  flatBase.name = stripext(flatBase.name);

  % install it
  disp(sprintf('(makeFlat) installing new flat base anatomy: %s', params.flatFileName));
  viewSet(view, 'newBase', flatBase);
  refreshMLRDisplay(viewNum);

  % remove the temporary off file (actually we should leave the
  % file since we need it to use the mrFlatViewer on the patch
  % later - but is it worth storing this some place in the
  % structure rather than as a file which migth get lost)? -jg
  %tempFileName = fullfile(params.path,params.flatFileName);
  %if isfile(tempFileName),delete(tempFileName);end
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%   myCutAndFlatten   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function flatSurface = myCutAndFlatten(params)

% set the name of the patch to cut
params.patchFileName = sprintf('%s_Patch.off',stripext(params.flatFileName));

% set lib path
mylibs = getenv('DYLD_LIBRARY_PATH');
setenv('DYLD_LIBRARY_PATH', sprintf('%s:/Users/eli/src/TFI/sw/lib/', mylibs));

% check for the SurfRelax program called 'surfcut'
[statusCut,result] = system('surfcut');
[statusFlat,result] = system('FlattenSurface.tcl');
if any([statusCut statusFlat]) ~= 0
  disp(sprintf('(makeFlat) Could not run the SurfRelax program surfcut. Make sure that you have SurfRelax correctly installed on your system.'));
  return;
else

  % cut and flatten
  degenFlag = 1; 
  distanceInc = 0;
  
  while degenFlag ~= 0
    disp(sprintf('(makeFlat) Cutting patch with radius of %i mm', params.radius+distanceInc))
    % cut the patch from the 3d mesh
    system(sprintf('surfcut -vertex %i -distance %i %s %s', ...
                   params.startVertex-1, params.radius+distanceInc, ... 
                   fullfile(params.path, params.innerCoordsFileName), ...
                   fullfile(params.path, params.patchFileName)));
    disppercent(inf);
    
    % flatten the patch
    disppercent(-inf, sprintf('(makeFlat) Flattening surface'));
    [degenFlag result] = system(sprintf('FlattenSurface.tcl %s %s %s', ...
                                        fullfile(params.path, params.outerCoordsFileName), ...
                                        fullfile(params.path, params.patchFileName), ...
                                        fullfile(params.path, params.flatFileName)));
    disppercent(inf);
    
    % if FlattenSurface failed, most likely b/c surfcut made a bad patch
    % increase the distance by one and try again.
    if (degenFlag ~= 0 ) && askuser('(makeFlat) Patch is degenerate, should I increase distance and reflatten?')
      distanceInc = distanceInc + 1;
    else
      degenFlag = 0;
      return;
    end

  end
  
end

% remove the 3d patch, b/c it isn't need for anything
if isfile(fullfile(params.path, params.patchFileName))
  system(sprintf('rm -rf %s', fullfile(params.path, params.patchFileName)));
end

% load the flat surface just created
flatSurface = importFlatOFF(params);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the surfaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [surf, params] = loadSurfHandler(params)
% we have already loaded the flat patch
% now load the rest of the surfaces

% read in the anatomy file
[surf.anat.data  surf.anat.hdr] = mlrImageReadNifti(fullfile(params.path, params.anatFileName));
% get vol2tal and vol2mag from the anatomy file
matFileName = [stripext(params.anatFileName) '.mat'];
if(exist([params.path '/' matFileName]))
  load([params.path '/' matFileName]);
  [tf base] = isbase(base);
  params.vol2mag = base.vol2mag;
  params.vol2tal = base.vol2tal;
  clear base
else
  params.vol2mag = [];
  params.vol2tal = [];
end

% load the white matter surface
surf.inner = loadSurfOFF(fullfile(params.path, params.innerCoordsFileName));
surf.inner = xformSurfaceWorld2Array(surf.inner, surf.anat.hdr);

% load the gray matter surface
surf.outer = loadSurfOFF(fullfile(params.path, params.outerCoordsFileName));
surf.outer = xformSurfaceWorld2Array(surf.outer, surf.anat.hdr);

% % read in the curvature file
[surf.curv, hdr] = loadVFF(fullfile(params.path, params.curvFileName));
% needs to be transposed to match the order of the vertices
surf.curv = surf.curv';           

% Extract permutation matrix to keep track of slice orientation.
surf.anat.permutationMatrix = getPermutationMatrix(surf.anat.hdr);

return

%%%%%%%%%%%%%%%%%%%%%%%
%%   runMrFlatMesh   %%
%%%%%%%%%%%%%%%%%%%%%%%
function[surf, params] = runMrFlatMesh(params)

% load the surfaces
[surf, params] = loadSurfHandler(params);

% project the surface out to an intermediate cortical deapth
corticalDepth = 0.5;
mesh.vertices = surf.inner.vtcs+corticalDepth*(surf.outer.vtcs-surf.inner.vtcs);

% create the mesh structure that mrFlatMesh expects
mesh.faceIndexList  = surf.inner.tris;
mesh.rgba           = surf.curv;
mesh.normal = surf.inner.vtcs - surf.outer.vtcs;

% run a modified version of the mrFlatMesh code
% this outputs and flattened surface
disppercent(-inf,'(makeFlat) Calling flattenSurfaceMFM');
surf.flat = flattenSurfaceMFM(mesh, [params.x params.y params.z], params.radius);
disppercent(inf);

% we need to figure out whether the flattened patch has been flipped
% during flattening
patch2parent = surf.flat.vertsToUnique(surf.flat.insideNodes);
vIn   = surf.inner.vtcs(patch2parent,:);
vOut  = surf.outer.vtcs(patch2parent,:);
vFlat = surf.flat.locs2d;
f     = surf.flat.uniqueFaceIndexList;

% this is the command to view the patch, in 2D or 3D
%hp = patch('vertices', v, 'faces', f, 'facecolor','none','edgecolor','black');

% loop through all of the faces
disppercent(-inf,'Checking winding direction');
wrapDir = zeros(1,length(f));wrapDirFlat = zeros(1,length(f));
for iFace = 1:length(f);
  % grab a triangle for inner 3D suface
  triIn = vIn(f(iFace,:),:);
  % grab a triangle for outer 3D suface
  triOut = vOut(f(iFace,:),:);
  % calculate the vector normal to the center of the two triangles
  triNorm = (mean(triIn) - mean(triOut)) + mean(triIn);
  % this is a formula that takes the vertices of the triangle and
  % computes the winding direction relative to a fourth point, which
  % is in this case the normal.  i.e. do the vertices of the triangle
  % go in a CW or CCW direction with respect to the normal. If this
  % determinant is positive the direction is CW and if the determinant
  % is negative it is CCW. see wikipedia:
  % http://en.wikipedia.org/wiki/Orientation_(topology)
  wrapDir(iFace) = det([cat(2,triIn, [1 1 1]'); triNorm 1]);
  
  % now do the same for the triangles in the flatpatch
  triFlat = vFlat(f(iFace,:),:);
  % in the flat patch, the z-dimension is always 0
  triFlat(:,3) = 0;
  % we want to compute the winding direction from above the surface
  % (i.e., the direction that we are viewing the surface from)
  triFlatNorm = [0 0 1];
  % same formula as above
  wrapDirFlat(iFace) = det([cat(2,triFlat, [1 1 1]'); triFlatNorm 1]);
  disppercent(iFace/length(f));
end
disppercent(inf);

% now check to see if the winding directions for the flat patch and 3D
% surface are the same or different. Note that because of the
% flattening process, some of the triangles may switch winding
% direction. But whether we should view the patch from above or below
% is decided by which view produces the least number of mismatches in
% winding direction.
match =    sum( sign(wrapDir) == sign(wrapDirFlat) );
misMatch = sum( sign(wrapDir) ~= sign(wrapDirFlat) );

if misMatch > match
  disp(sprintf('(makeFlat) The patch is NOT oriented properly (misMatch=%0.2f%%), X-Y flipping the patch...',100*misMatch/(match+misMatch)));
  surf.flat.locs2d = cat(2, surf.flat.locs2d(:,2), surf.flat.locs2d(:,1));
else
  disp(sprintf('(makeFlat) The patch is oriented properly (misMatch=%0.2f%%), not going to flip it...',100*misMatch/(match+misMatch)));
end

% save the patch
writePatchOFF(surf, params);

% convert surface into a flat patch base strucutre, by loading it back up
surf = importFlatOFF(params);

return;



% writeOFF.m
%
%      usage: writeOFF()
%         by: eli merriam
%       date: 10/25/07
%    purpose: 
%
function retval = writePatchOFF(surf, params)

% check arguments
if ~any(nargin == [ 0 1 2])
  help writeOFF
  return
end

% Vertices
vertices = [surf.flat.locs2d(:,1) surf.flat.locs2d(:,2)];
vertices = cat(1, vertices', zeros(1, length(vertices)));

% triangles(1) is number of vert/triangle: 3
% triangles(2:4) are the vertices of the triangles
% triangles(5) is color: 0
triangles =  surf.flat.uniqueFaceIndexList'-1;
triangles = cat(1, ones(1,length(triangles))*3, triangles, zeros(1,length(triangles),1));

patch2parent = surf.flat.vertsToUnique(surf.flat.insideNodes);

% write the OFF format file 
fid = fopen(fullfile(params.path, params.flatFileName), 'w', 'ieee-be');
if fid == -1
  disp(sprintf('(makeFlat) Could not open %s for saving the flat patch. Do you have permissions set correctly in the directory?',fullfile(params.path, params.flatFileName)));
  disp(sprintf('(makeFlat) Stopping execution. Type dbquit to get out of keyboard mode'));
  keyboard
end
fprintf(fid, '#PATCH\n');
fprintf(fid, '#parent_surface=%s\n', fullfile(params.path, params.innerCoordsFileName));
fprintf(fid, '#parent_dimensions=%i %i %i\n', surf.inner.Nvtcs,  surf.inner.Ntris, 1);
fprintf(fid, '#patch_dimensions=%i %i %i\n', length(vertices), length(triangles), 1 );
fprintf(fid, '#parent_vertex_indexes:\n');
for i=1:length(surf.flat.uniqueVertices)
  fprintf(fid, '#%i %i\n', i-1, patch2parent(i)-1);
end

fprintf(fid, 'OFF BINARY\n');
fwrite(fid, [size(vertices,2) size(triangles,2) 0], 'int32'); 

% Vertices
fwrite(fid, vertices, 'float32');

% Faces
fwrite(fid, triangles, 'int32');

% Close file
fclose(fid);

return;
