% importSurfaceOFF.m
%
%      usage: v = importSurfaceOFF
%         by: justin gardner
%       date: 10/24/07
%    purpose: Import a pair of inner and outer cortical surfaces in OFF format
%
%    base = importSurfaceOFF('/path/to/surfaces/subject_left_WM.off', bothHemiFlag);
%    viewSet(getMLRView, 'newbase', base);
%
function base = importSurfaceOFF(pathStr,bothHemiFlag)

% check arguments
if ~any(nargin == [0 1 2])
  help importSurfaceOFF
  return
end
base = [];
if ieNotDefined('bothHemiFlag'), 
  bothHemiFlag = 0;
  disp(sprintf('Only loading surfaces for one hemisphere'))
else
  disp(sprintf('Loading surfaces for both hemispheres'))
end

if ieNotDefined('pathStr')
  % Open dialog box to have user choose the file
  startPathStr = mrGetPref('volumeDirectory');
  filterspec = {'*.off','Off Surface file (*.off)'};
  title = 'Choose outer surface file';
  pathStr = mlrGetPathStrDialog(startPathStr,title,filterspec,'off');
end

% Aborted
if ieNotDefined('pathStr'),return,end

% get surface name using mrSurfViewer
[filepath filename] = fileparts(pathStr);
thispwd = pwd;
if ~isempty(filepath),cd(filepath);end
params1 = mrSurfViewer(filename);

% Aborted
if isempty(params1),cd(thispwd);return,end

% Create the base
base.hdr = mlrImageReadNiftiHeader(params1.anatomy);
if isempty(base.hdr)
  mrWarnDlg(sprintf('(imortSurfaceOFF) Could not load anatomy file: %s',params1.anatomy));
  base = [];
  return
end
base.name = filename;
base.permutationMatrix = getPermutationMatrix(base.hdr);
base.range = [-1.5 1.5];
base.clip = [-1.5 1.5];

% get the vol2mag and vol2tal fields from the volume anatomy
% do it using a subfunction so don't confuse the two 'base' structure variables
anatomyMatFile = sprintf('%s.mat',stripext(params1.anatomy));
base.vol2tal = getBaseField(anatomyMatFile,'vol2tal');
base.vol2mag = getBaseField(anatomyMatFile,'vol2mag');

% load both hemispheres
if bothHemiFlag
  % get the params for the other hemisphere
  % was a left hemisphere passed?
  leftFlag = strfind(lower(params1.innerSurface), 'left');
  if leftFlag
    prefix = params1.innerSurface(1:leftFlag-1);
    postfix = params1.innerSurface(leftFlag+4:end);
    rightFilename = sprintf('%sright%s', prefix, postfix);
    params2 = mrSurfViewer(rightFilename);  
  end
  % or was it a right hemispehre
  rightFlag = strfind(lower(params1.innerSurface), 'right');
  if rightFlag
    prefix = params1.innerSurface(1:rightFlag-1);
    postfix = params1.innerSurface(rightFlag+5:end);
    leftFilename = sprintf('%sleft%s', prefix, postfix);
    params2 = mrSurfViewer(leftFilename);  
  end
  % and add the second curvature
  data1(1,:,1) = loadVFF(params1.curv);
  data2(1,:,1) = loadVFF(params2.curv);
  base.data(1,:,1) = cat(2, data1, data2);
  % load both inner surfaces
  innerSurface1 = loadSurfOFF(params1.innerSurface);
  innerSurface2 = loadSurfOFF(params2.innerSurface);
  innerSurface = combineSurfaces(innerSurface1, innerSurface2);
  % load both outer surfaces
  outerSurface1 = loadSurfOFF(params1.outerSurface);
  outerSurface2 = loadSurfOFF(params2.outerSurface);
  outerSurface = combineSurfaces(outerSurface1, outerSurface2);
  % load the inner coords   
  if strcmp(params1.innerCoords,'Same as surface')
    inner = innerSurface;
  else
    inner1 = loadSurfOFF(params1.innerCoords);
    inner2 = loadSurfOFF(params2.innerCoords);
    inner = combineSurfaces(inner1,inner2);
  end
  if strcmp(params1.outerCoords,'Same as surface')
    outer = outerSurface;
  else
    outer1 = loadSurfOFF(params1.outerCoords);
    outer2 = loadSurfOFF(params2.outerCoords);
    outer = combineSurfaces(inner1,inner2);
  end

else                                    
  % or else a single hemisphere...
  % load the curvature
  base.data(1,:,1) = loadVFF(params1.curv);
  % load inner surfaces
  innerSurface = loadSurfOFF(params1.innerSurface);
  % load outer surface
  outerSurface = loadSurfOFF(params1.outerSurface);
  if strcmp(params1.innerCoords,'Same as surface')
    inner = innerSurface;
  else
    inner = loadSurfOFF(params1.innerCoords);
  end
  if strcmp(params1.outerCoords,'Same as surface')
    outer = outerSurface;
  else
    outer = loadSurfOFF(params1.outerCoords);
  end
end  

% now convert all surfaces to array coordinates
innerSurface = xformSurfaceWorld2Array(innerSurface,base.hdr);
outerSurface = xformSurfaceWorld2Array(outerSurface,base.hdr);
inner = xformSurfaceWorld2Array(inner,base.hdr);
outer = xformSurfaceWorld2Array(outer,base.hdr);

% save names
base.coordMap.path = filepath;
base.coordMap.innerSurfaceFileName = params1.innerSurface;
base.coordMap.innerCoordsFileName = params1.innerCoords;
base.coordMap.outerSurfaceFileName = params1.outerSurface;
base.coordMap.outerCoordsFileName = params1.outerCoords;
base.coordMap.curvFileName = params1.curv;
base.coordMap.anatFileName = params1.anatomy;
% save coords
base.coordMap.innerCoords(1,:,1,1)  = inner.vtcs(:,1);
base.coordMap.innerCoords(1,:,1,2)  = inner.vtcs(:,2);
base.coordMap.innerCoords(1,:,1,3)  = inner.vtcs(:,3);
base.coordMap.innerVtcs = innerSurface.vtcs;
base.coordMap.tris = innerSurface.tris;
base.coordMap.outer = params1.outerSurface;
base.coordMap.outerCoords(1,:,1,1)  = outer.vtcs(:,1);
base.coordMap.outerCoords(1,:,1,2)  = outer.vtcs(:,2);
base.coordMap.outerCoords(1,:,1,3)  = outer.vtcs(:,3);
base.coordMap.outerVtcs = outerSurface.vtcs;
base.coordMap.coords = base.coordMap.innerCoords;
base.coordMap.dims = base.hdr.dim([3 2 4])';
base.type = 2;

cd(thispwd);

function val = getBaseField(matFilename,fieldname)
if ~exist(matFilename,'file') % if the anatomy doesn't have these fields set
  val = [];
else % if the volume anatomy has the fields, use them
  load(matFilename);
  eval(sprintf('val = base.%s;',fieldname))
end


function surf = combineSurfaces(surf1, surf2);
surf.filename = sprintf('%_%s', surf1.filename, surf2.filename);
surf.Nvtcs = surf1.Nvtcs + surf2.Nvtcs;
surf.Ntris = surf1.Ntris + surf2.Ntris;
surf.Nedges = 0;
surf.vtcs = [surf1.vtcs; surf2.vtcs];
surf.tris = [surf1.tris; surf2.tris+surf1.Nvtcs];
