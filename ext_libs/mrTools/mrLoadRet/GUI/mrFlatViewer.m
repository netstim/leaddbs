% mrFlatViewer.m
%
%       $Id$	
%      usage: params = mrFlatViewer(flat,<outer>,<inner>,<curv>,<anat>,<viewNum>)
%         by: justin gardner, originally based on surfViewer by eli merriam
%       date: 10/09/07
%    purpose: Displays a flattened patch. Flat can be a name of a file:
%
%             mrFlatViewer('jg_left_MT_flat');
%
%             or can be a structure that specifies the location of a patch
%             defined by a point/radius:
%
%             flat.path = pwd;
%             flat.parentSurfaceName = 'jg_left_WM';
%             flat.startPoint = [200 50 100];
%             flat.radius = 50;
%             params = mrFlatViewer(flat);
%
function retval = mrFlatViewer(flat,outer,inner,curv,anat,viewNum)

% check arguments
if ~any(nargin == [1 2 3 4 5 6])
  help mrFlatViewer
  return
end
if nargout == 1
  retval = [];
end

% if passed in a string check to see if
% it needs an extension
if isstr(flat)
  if isfile(sprintf('%s.off',stripext(flat)));
    flat = sprintf('%s.off',stripext(flat));
  end
end

% see how we are being called
if (nargin == 1) && isstr(flat) && ~isfile(flat)
  event = flat;
else
  event = 'init';
  % set defaults
  if ieNotDefined('outer'),outer = {};end
  if ieNotDefined('inner'),inner = {};end
  if ieNotDefined('curv'),curv = {};end
  if ieNotDefined('anat'),anat = {};end
  if ieNotDefined('viewNum'),viewNum = [];end
  % make everybody a cell array
  flat = cellArray(flat);
  outer = cellArray(outer);
  inner = cellArray(inner);
  curv = cellArray(curv);
  anat = cellArray(anat);
end

switch (event)
 case 'init'
  retval = initHandler(flat,outer,inner,curv,anat,viewNum);
 case {'vSlider','hSlider'}
  sliderHandler;
 case {'edit'}
  editHandler;
 otherwise
  disp(sprintf('(mrFlatViewer) Could not find flat file %s',flat{1}));
end

%%%%%%%%%%%%%%%%%%%%%%
%%   init handler   %%
%%%%%%%%%%%%%%%%%%%%%%
function retval = initHandler(flat,outer,inner,curv,anat,viewNum)

global gFlatViewer;
gFlatViewer = [];
gFlatViewer.mismatchWarning = 0;
retval = [];
disppercent(-inf,'(mrFlatViewer) Loading surfaces');

% load the flat
if isstr(flat{1})
  [flatPath flat{1}] = fileparts(sprintf('%s.off',stripext(flat{1})));
  flatdir = dir(fullfile(flatPath,'*.off'));
  gFlatViewer.path = flatPath;

  gFlatViewer.flat = loadSurfOFF(fullfile(flatPath, sprintf('%s.off', stripext(flat{1}))));
  if isempty(gFlatViewer.flat) || ~isfield(gFlatViewer.flat,'parentSurfaceName');
    disp(sprintf('(mrFlatViewer) %s is not a flat file',flat{1}));
    return
  end
  % remove any paths
  gFlatViewer.flat.parentSurfaceName = getLastDir(gFlatViewer.flat.parentSurfaceName);

elseif isfield(flat{1},'radius')
  % if this is a structure, then we are being called from makeFlat
  % with coordinates and a radius
  flatdir = [];
  gFlatViewer.path = flat{1}.path;
  flatPath = flat{1}.path;
  %first convert coordinates from current base to the surface base 
  if isempty(anat)
    [filename, pathname] = uigetfile({'*.hdr;*.nii','Nifti file (*.hdr/*.nii)'},'Select 3D Anatomy File',flat{1}.path);
    anat{1} = [pathname filename];
    anatomyFile=anat{1};
  else
  % (we assume the base anatomy is in the same folder as the parent surface)
    anatomyFile=fullfile(flat{1}.path,anat{1});
   end
  hdr= cbiReadNiftiHeader(anatomyFile);
  baseStartPoint = hdr.sform44 \ viewGet(viewNum,'basexform') * [flat{1}.startPoint';1];
  gFlatViewer.flat = makeFlatFromRadius(flat{1},flat{1}.radius,baseStartPoint(1:3)',flat{1}.parentSurfaceName);
  if isempty(gFlatViewer.flat),return,end
  flat{1} = gFlatViewer.flat.name;
elseif isfield(flat{1},'vtcs')
  % is a passed in flat 
  gFlatViewer.flat = flat{1};
  flat{1} = gFlatViewer.flat.filename;
  flatdir = [];
  gFlatViewer.path = gFlatViewer.flat.path;
  flatPath = gFlatViewer.path;
end

% look for flats with same parent
for i = 1:length(flatdir)
  if ~isempty(strcmp(lower(flatdir(i).name),'flat')) || ~isempty(strcmp(lower(flatdir(i).name),'patch'))
    flatfile = loadSurfOFF(fullfile(flatPath, flatdir(i).name),1);
    if isfield(flatfile,'parentSurfaceName')
      if strcmp(flatfile.parentSurfaceName,gFlatViewer.flat.parentSurfaceName)
	if ~strcmp(flatdir(i).name,flat)
	  flat{end+1} = flatdir(i).name;
	end
      end
    end
  end
end

% load up the surfaces
checkForMore = 1;
if isempty(inner)
  % guess the names
  if isfile(fullfile(flatPath,gFlatViewer.flat.parentSurfaceName))
    inner{1} = gFlatViewer.flat.parentSurfaceName;
  else
    % go look for it
    [filename pathname] = uigetfile({'*.off','Surface file (*.off)'},sprintf('Find parent surface %s',gFlatViewer.flat.parentSurfaceName),flatPath);
    if isempty(filename),return,end
    filename = fullfile(getRelativePath(flatPath,pathname),filename);
    % save it
    gFlatViewer.flat.parentSurfaceName = filename;
    % set mismatch warning on, so we check to see if there is a mismatch
    gFlatViewer.mismatchWarning = 1;
    checkForMore = 0;
    inner{1} = filename;
  end
end

% guess anything with the right stem
if checkForMore
  innerDir = dir(sprintf('%s*.off',fullfile(flatPath,stripext(stripext(inner{1}),'WM'))));
  for i = 1:length(innerDir)
    % don't choose anything we already have or with GM, flat or patch in the title
    if ~any(strcmp(innerDir(i).name,inner)) && isempty(strfind(innerDir(i).name,'GM')) && (isempty(strfind(lower(innerDir(i).name),'flat')) || ~isempty(strfind(lower(innerDir(i).name),'inflate'))) && isempty(strfind(lower(innerDir(i).name),'patch'))
      % check to make sure it has the correct number of vertices
      surf = myLoadSurface(innerDir(i).name, flatPath,1);
      if ~isempty(surf)
	inner{end+1} = innerDir(i).name;
      end
    end
  end
end

% now try to find the first loadable one
for i = 1:length(inner)
  gFlatViewer.surfaces.inner = myLoadSurface(fullfile(flatPath, inner{i}));
  gFlatViewer.mismatchWarning = 0;
  if ~isempty(gFlatViewer.surfaces.inner),break,end
end
% if we didn't load anything then quit
if isempty(gFlatViewer.surfaces.inner)
  mrWarnDlg(sprintf('(mrFlatViewer) Could not load inner surface %s',inner{1}));
  return
else
  inner = putOnTopOfList(inner{i},inner);
end
inner{end+1} = 'Find file';

% load the outer surface
if isempty(outer)
  % if we weren't passed in anything try to find them
  filename = sprintf('%sGM.off',stripext(stripext(inner{1}),'WM'));
  if isfile(fullfile(flatPath,filename))
    outer{1} = filename;
  else
    % go look for it
    [filename pathname] = uigetfile({'*.off','Surface file (*.off)'},sprintf('Find outer surface'),flatPath);
    if isempty(filename),return,end
    filename = fullfile(getRelativePath(flatPath,pathname),filename);
    % set mismatch warning on, so we check to see if there is a mismatch
    gFlatViewer.mismatchWarning = 1;
    checkForMore = 0;
    outer{1} = filename;
  end
end
% guess anything with the right stem
if checkForMore
  outerDir = dir(sprintf('%s*.off',fullfile(flatPath,stripext(stripext(inner{1}),'WM'))));
  for i = 1:length(outerDir)
    % don't choose anything we already have or with WM, flat or patch in the title
    if ~any(strcmp(innerDir(i).name,inner)) && isempty(strfind(innerDir(i).name,'WM')) && (isempty(strfind(lower(innerDir(i).name),'flat')) || ~isempty(strfind(lower(innerDir(i).name),'inflate'))) && isempty(strfind(lower(innerDir(i).name),'patch'))
      % check to make sure it has the correct number of vertices
      surf = myLoadSurface(outerDir(i).name, flatPath, 1);
      if ~isempty(surf)
	outer{end+1} = outerDir(i).name;
      end
    end
  end
end

% now try to find the first loadable one
for i = 1:length(outer)
  gFlatViewer.surfaces.outer = myLoadSurface(sprintf('%s',outer{i}), flatPath);
  if ~isempty(gFlatViewer.surfaces.outer),break,end
end
% if we didn't load anything then quit
if isempty(gFlatViewer.surfaces.outer)
  return
else
  outer = putOnTopOfList(outer{i},outer);
end
outer{end+1} = 'Find file';

% load the curvature
checkForMore = 1;
if isempty(curv)
  curvGuess = sprintf('%s_Curv.vff', fullfile(flatPath, stripext(inner{1})));
  secondCurvGuess = sprintf('%sCurv.vff', fullfile(flatPath, stripext(stripext(inner{1}),'WM')));
  if isfile(curvGuess)
    curv{1} = getLastDir(curvGuess);
  elseif isfile(secondCurvGuess)
    curv{1} = getLastDir(secondCurvGuess);
  else
    % go look for it
    [filename pathname] = uigetfile({'*.vff','Curvature file (*.vff)'},sprintf('Find curvature file'),flatPath);
    if isempty(filename),return,end
    filename = fullfile(getRelativePath(flatPath,pathname),filename);
    % set mismatch warning on, so we check to see if there is a mismatch
    gFlatViewer.mismatchWarning = 1;
    checkForMore = 0;
    curv{1} = filename;
  end 
end
if checkForMore
  % add any vff file
  curvDir = dir(sprintf('%s/*.vff', flatPath));
  for i = 1:length(curvDir)
    if ~any(strcmp(curvDir(i).name,curv))
      % check length of file matches our patch
      vffhdr = myLoadCurvature(curvDir(i).name, flatPath, 1);
      if ~isempty(vffhdr)
	curv{end+1} = curvDir(i).name;
      end
    end
  end
end

for i = 1:length(curv)
  gFlatViewer.curv = myLoadCurvature(sprintf('%s', fullfile(flatPath, curv{1})));
  if ~isempty(gFlatViewer.curv),break,end
end
% if we didn't load anything then quit
if isempty(gFlatViewer.curv)
  mrWarnDlg(sprintf('(mrFlatViewer) Could not find a matching .vff curvature file. You will need to use calcCurvature to compute a curvature file from the inner (WM) and outer (GM) surfaces.',outer{1},stripext(outer{1})));
  return
else
  curv = putOnTopOfList(curv{i},curv);
end
disppercent(inf);

% from now on, complain for mismatch of surface nodes and patches
gFlatViewer.mismatchWarning = 1;

% guess any nifti file for anatomy
anatDir = [dir(sprintf('%s/*.hdr', flatPath)); dir(sprintf('%s/*.nii', flatPath))];
for i = 1:length(anatDir)
  if ~any(strcmp(anatDir(i).name,anat))
    anat{end+1} = anatDir(i).name;
  end
end

% check for 'canonical hdr'
anatCanonicalDir = [dir(sprintf('%s/../*.hdr', flatPath)); dir(sprintf('%s/../*.nii', flatPath))];
for i= 1:length(anatCanonicalDir)
  anat = putOnTopOfList(fullfile('..',anatCanonicalDir(i).name),anat);
end

% now try to open the file
% if we haven't yet found a good candidate, then ask user
if isempty(anat)
  % go look for it
  [filename pathname] = uigetfile({'*.hdr;*.nii','Nifti file (*.hdr/*.nii)'},'Find 3D anatomy',flatPath);
  if isnumeric(filename),return,end
  filename = fullfile(getRelativePath(flatPath,pathname),filename);
  anat{1} = filename;
end

% now go through and load (to make sure we have valid nifti files
validAnat = [];
for iAnat = 1:length(anat)
  thisAnatName = fullfile(flatPath, anat{iAnat});
  % assume not valid at first
  validAnat(iAnat) = false;
  % check for file
  if isfile(thisAnatName)
    [gFlatViewer.anat.data gFlatViewer.anat.hdr] = mlrImageReadNifti(thisAnatName);
    % if no data, then there was a failure to load
    if isempty(gFlatViewer.anat.data)
      mrWarnDlg(sprintf('(mrFlatViewer) Could not load nifti file %s',thisAnatName));
      validAnat(iAnat) = false;
    else
      gFlatViewer = xformSurfaces(gFlatViewer);
      validAnat(iAnat) = true;
    end
  end
end

% remove all non-working anatomies (also flip order, so that the last one loaded above
% is the one we put at the top of the list - so that the list matches what is loaded).
anat = {anat{fliplr(find(validAnat))}};
if isempty(anat)      
  % there are no candidate anatomies. This is a problem.
  mrWarnDlg('(mrFlatViewer) Could not find any valid 3D anatomies');
  return
end

% save the view
gFlatViewer.viewNum = viewNum;

% select the window
gFlatViewer.f = selectGraphWin;
set(gFlatViewer.f,'renderer','OpenGL');

% positions on figure
figLeft = 10;figBottom = 10;
sliderWidth = 20;sliderLength = 200;spacer = 10;
editWidth = 40;editHeight = 20;

% set up horizontal and vertical slider
gFlatViewer.hSliders.v = uicontrol('Style','slider','Position',[figLeft figBottom+sliderWidth sliderWidth sliderLength],'Min',-180,'Max',180,'SliderStep',[15 45]./360,'Callback','mrFlatViewer(''vSlider'')','TooltipString','Rotate around y-axis');
gFlatViewer.hSliders.vText = uicontrol('Style','Edit','Position',[figLeft figBottom+sliderWidth+sliderLength+spacer editWidth editHeight],'Callback','mrFlatViewer(''edit'')','String','0','HorizontalAlignment','Center');
gFlatViewer.hSliders.h = uicontrol('Style','slider','Position',[figLeft+sliderWidth figBottom sliderLength sliderWidth],'Min',-180,'Max',180,'SliderStep',[15 45]./360,'Callback','mrFlatViewer(''hSlider'')','TooltipString','Rotate around z-axis');
gFlatViewer.hSliders.hText = uicontrol('Style','Edit','Position',[figLeft+sliderLength+3*spacer figBottom editWidth editHeight],'Callback','mrFlatViewer(''edit'')','String','0');


% set they we are viewing white matter
gFlatViewer.whichSurface = 1;
gFlatViewer.patchColoring = 1;
gFlatViewer.displayROIs = 0;
% and display surface
dispSurface;
setViewAngle(0,0);

editable = 0;

% set up the parameters
paramsInfo = {};
gFlatViewer.guiloc.whichSurface = 1;
gFlatViewer.guiloc.filenames = 3;
% Now give choice of viewing gray or white
gFlatViewer.whichSurfaceTypes = {'Outer (Gray matter) surface','Inner (White matter) surface','3D Anatomy','Patch'};
paramsInfo{end+1} = {'whichSurface',gFlatViewer.whichSurfaceTypes,'type=popupmenu','callback',@whichSurfaceCallback,'Choose which surface to view the patch on'};
gFlatViewer.patchColoringTypes = {'Uniform','Right in red','Rostral in red','Dorsal in red','Positive curvature in red','Negative curvature in red','Compressed areas in red','Stretched areas in red','High outer areal distortion in red','High inner areal distortion in red'};
if ~isempty(gFlatViewer.viewNum)
  gFlatViewer.patchColoringTypes{end+1} = 'Current overlay';
  gFlatViewer.patchColoringTypes{end+1} = 'Current overlay with patch';
end
gFlatViewer.patchColoringTypes{end+1} = 'None';
paramsInfo{end+1} = {'patchColoring',gFlatViewer.patchColoringTypes,'type=popupmenu','Choose how to color the patch','callback',@patchColoringCallback};
if ~isempty(gFlatViewer.viewNum)
  gFlatViewer.guiloc.filenames = gFlatViewer.guiloc.filenames+1;
  paramsInfo{end+1} = {'displayROIs',0,'type=checkbox','Display the ROIs','callback',@whichSurfaceCallback};
end
paramsInfo{end+1} = {'path', flatPath,'editable=0','The directory path to the flat file'};
if isfield(gFlatViewer.flat,'radius')
  paramsInfo{end+1} = {'flatFileName',flat{1},'editable=1','The flat patch file'};
  paramsInfo{end+1} = {'radius',gFlatViewer.flat.radius,'incdec=[-5 5]','minmax=[1 inf]','callback',@setFlatRadius,'Set the radius in mm of the flat patch'};
  paramsInfo{end+1} = {'x',gFlatViewer.flat.startPoint(1),'incdec=[-10 10]','minmax=[1 inf]','callback',@setFlatStartPoint,'Set the start x position of patch. Note that if you modify this field it will get reset to the closest [x y s] point that is on the surface.'};
  paramsInfo{end+1} = {'y',gFlatViewer.flat.startPoint(2),'incdec=[-10 10]','minmax=[1 inf]','callback',@setFlatStartPoint,'Set the start y position of patch. Note that if you modify this field it will get reset to the closest [x y s] point that is on the surface.'};
  paramsInfo{end+1} = {'z',gFlatViewer.flat.startPoint(3),'incdec=[-10 10]','minmax=[1 inf]','callback',@setFlatStartPoint,'Set the start z position of patch. Note that if you modify this field it will get reset to the closest [x y s] point that is on the surface.'};
elseif ~editable && (length(flat) == 1)
  paramsInfo{end+1} = {'flatFileName',flat{1},'editable=0','The flat patch file'};
else
  paramsInfo{end+1} = {'flatFileName',flat,'The flat patch file','callback',@switchFlat};
end
if ~editable && (length(outer) == 1)
  paramsInfo{end+1} = {'outerCoordsFileName',outer{1},'editable=0','The outer (gray matter) file'};
else
  paramsInfo{end+1} = {'outerCoordsFileName',outer,'The outer (gray matter) file','callback',@switchFile,'callbackArg=outerCoordsFileName'};
end
if ~editable && (length(inner) == 1)
  paramsInfo{end+1} = {'innerCoordsFileName',inner{1},'editable=0','The inner (white matter) file'};
else
  paramsInfo{end+1} = {'innerCoordsFileName',inner,'The inner (white matter) file','callback',@switchFile,'callbackArg=innerCoordsFileName'};
end
if ~editable && (length(curv) == 1)
  paramsInfo{end+1} = {'curvFileName',curv{1},'editable=0','The curvature file. This is a file that can be created from the inner (WM) and outer (GM)  surface with the command calcCurvature'};
else
  paramsInfo{end+1} = {'curvFileName',curv,'The curvature file','callback',@switchFile,'callbackArg=curvFileName'};
end
if ~editable && (length(anat) == 1)
  paramsInfo{end+1} = {'anatFileName',anat{1},'editable=0','The 3D anatomy file'};
else
  paramsInfo{end+1} = {'anatFileName',anat,'The 3D anatomy file','callback',@switchAnatomy};
end

% put up dialog
if isfield(gFlatViewer.flat,'radius')
  params = mrParamsDialog(paramsInfo,'Set parameters for flat patch');
else
  params = mrParamsDialog(paramsInfo,'View flat patch location on surface');
end
if isempty(params)
  retval = [];
else
  params.flatFileName = setext(fixBadChars(stripext(params.flatFileName)),'off');
  % return also the startVertex if this is a radius /start position
  if isfield(gFlatViewer.flat,'startVertex')
    params.startVertex = gFlatViewer.flat.startVertex;
  end
  retval = params;
end
if ishandle(gFlatViewer.f)
  close(gFlatViewer.f);
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   patchColoringCallback   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function patchColoringCallback(params)

global gFlatViewer
gFlatViewer.patchColoring = find(strcmp(params.patchColoring,gFlatViewer.patchColoringTypes));
if gFlatViewer.whichSurface == 3
  hPos = round(get(gFlatViewer.hSliders.h,'Value'));
  dispVolume(3,hPos);
else
  dispSurface;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   whichSurfaceCallback   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function whichSurfaceCallback(params)

global gFlatViewer;

% set the roi drawing
if ~isfield(params,'displayROIs')
  params.displayROIs = gFlatViewer.displayROIs;
end

% get which surface to draw
lastWhichSurface = gFlatViewer.whichSurface;

refreshFlatViewer(find(strcmp(params.whichSurface,gFlatViewer.whichSurfaceTypes)),params.displayROIs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   refreshFlatViewer   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function refreshFlatViewer(whichSurface,displayROIs,force)

global gFlatViewer;

if ieNotDefined('whichSurface')
  whichSurface = gFlatViewer.whichSurface;
end
if ieNotDefined('displayROIs')
  displayROIs = gFlatViewer.displayROIs;
end
if ieNotDefined('force')
  force = 0;
end
% what surface/rois are being displayed now
lastWhichSurface = gFlatViewer.whichSurface;
lastDisplayROIs = gFlatViewer.displayROIs;

% get which surface to draw
if force || (whichSurface ~= lastWhichSurface) || (lastDisplayROIs ~= displayROIs)
  % set which surface and display
  gFlatViewer.whichSurface = whichSurface;
  gFlatViewer.displayROIs = displayROIs;  
  % 1,2 are surfaces
  if whichSurface <= 2
    % if we are displaying the 3D anatomy, 
    % then switch to the surface view
    if lastWhichSurface > 2
      switchToSurface;
    else
      dispSurface;
    end
  % 3 is the volume
  elseif whichSurface == 3
    % switch to the volume view
    switchToVolume;
  else
    % switch to the volume view
    switchToFlat;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%   switchToSurface   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function switchToSurface

global gFlatViewer;
set(gFlatViewer.hSliders.v,'Visible','on');
set(gFlatViewer.hSliders.vText,'Visible','on');
set(gFlatViewer.hSliders.h,'Visible','on');
set(gFlatViewer.hSliders.hText,'Visible','on');
set(gFlatViewer.hSliders.h,'SliderStep',[15 45]./360);
set(gFlatViewer.hSliders.h,'Value',0);
set(gFlatViewer.hSliders.h,'Min',-180);
set(gFlatViewer.hSliders.h,'Max',180);
set(gFlatViewer.hSliders.h,'TooltipString','Rotate around z-axis');
set(gFlatViewer.hSliders.v,'Value',0);
set(gFlatViewer.hSliders.vText,'String',0);
set(gFlatViewer.hSliders.hText,'String',0);
dispSurface;
setViewAngle(0,0);

%%%%%%%%%%%%%%%%%%%%%%%%
%%   switchToVolume   %%
%%%%%%%%%%%%%%%%%%%%%%%%
function switchToVolume

global gFlatViewer;
initSlice = 127;
set(gFlatViewer.hSliders.h,'Visible','on');
set(gFlatViewer.hSliders.hText,'Visible','on');
set(gFlatViewer.hSliders.v,'Visible','off');
set(gFlatViewer.hSliders.vText,'Visible','off');
set(gFlatViewer.hSliders.h,'SliderStep',[1 16]./256);
set(gFlatViewer.hSliders.h,'Value',initSlice);
set(gFlatViewer.hSliders.h,'Min',1);
set(gFlatViewer.hSliders.h,'Max',256);
set(gFlatViewer.hSliders.h,'TooltipString','Change viewing slice');
set(gFlatViewer.hSliders.hText,'String',num2str(initSlice));
dispVolume(3,initSlice);

%%%%%%%%%%%%%%%%%%%%%%
%%   switchToFlat   %%
%%%%%%%%%%%%%%%%%%%%%%
function switchToFlat

global gFlatViewer;
initSlice = 127;
% turn off sliders if this is a real flattened flat
if ~isfield(gFlatViewer.flat,'radius')
  set(gFlatViewer.hSliders.v,'Visible','off');
  set(gFlatViewer.hSliders.vText,'Visible','off');
  set(gFlatViewer.hSliders.h,'Visible','off');
  set(gFlatViewer.hSliders.hText,'Visible','off');
  dispSurface;
else
  set(gFlatViewer.hSliders.v,'Visible','on');
  set(gFlatViewer.hSliders.vText,'Visible','on');
  set(gFlatViewer.hSliders.h,'Visible','on');
  set(gFlatViewer.hSliders.hText,'Visible','on');
  set(gFlatViewer.hSliders.h,'SliderStep',[15 45]./360);
  set(gFlatViewer.hSliders.h,'Value',0);
  set(gFlatViewer.hSliders.h,'Min',-180);
  set(gFlatViewer.hSliders.h,'Max',180);
  set(gFlatViewer.hSliders.h,'TooltipString','Rotate around z-axis');
  set(gFlatViewer.hSliders.v,'Value',0);
  set(gFlatViewer.hSliders.vText,'String',0);
  set(gFlatViewer.hSliders.hText,'String',0);
  dispSurface;
  setViewAngle(0,0);
end

%%%%%%%%%%%%%%%%%%%%%%%
%%   sliderHandler   %%
%%%%%%%%%%%%%%%%%%%%%%%
function sliderHandler

global gFlatViewer;

% get slider position
hPos = round(get(gFlatViewer.hSliders.h,'Value'));
vPos = round(get(gFlatViewer.hSliders.v,'Value'));

% set the edit fields
set(gFlatViewer.hSliders.hText,'String',num2str(hPos));
set(gFlatViewer.hSliders.vText,'String',num2str(vPos));

if any(gFlatViewer.whichSurface == [1 2 4])
  setViewAngle(hPos,vPos);
else
  dispVolume(3,hPos);
end

%%%%%%%%%%%%%%%%%%%%%%%
%%   editHandler   %%
%%%%%%%%%%%%%%%%%%%%%%%
function editHandler

global gFlatViewer;

hPos = str2num(get(gFlatViewer.hSliders.hText,'String'));
vPos = str2num(get(gFlatViewer.hSliders.vText,'String'));

% make it fit into -180:180
hPos = round(mod(hPos+180,360)-180);
vPos = round(mod(vPos+180,360)-180);

% set slider position
set(gFlatViewer.hSliders.h,'Value',hPos);
set(gFlatViewer.hSliders.v,'Value',vPos);

% set the edit fields
set(gFlatViewer.hSliders.hText,'String',num2str(hPos));
set(gFlatViewer.hSliders.vText,'String',num2str(vPos));

if gFlatViewer.whichSurface <= 2
  setViewAngle(hPos,vPos);
else
  dispVolume(3,hPos);
end

%%%%%%%%%%%%%%%%%%%%%%
%%   setViewAngle   %%
%%%%%%%%%%%%%%%%%%%%%%
function setViewAngle(hPos,vPos)

% flip the sign to make rotations go in the "right" direction
hPos = -hPos;vPos = -vPos;
% somehow 90 and 180 are a problem for matlab
if abs(vPos) == 90,vPos = sign(vPos)*91;,end
if abs(hPos) == 90,hPos = sign(hPos)*91;,end
if abs(hPos) == 179,hPos = sign(hPos)*179;,end

% set the camera taret to center
camtarget([0 0 0]);

% set the size of the field of view in degrees
% i.e. 90 would be very wide and 1 would be ver
% narrow. 7 seems to fit the whole brain nicely
camva(7);

% set the view angle
view(hPos,vPos);

% change the camera position to avoid the volume
% flipping back and forth, another starnge matlab thing
if (vPos >= 90) || (vPos < -90)
  camup([0 0 -1]);
else
  camup([0 0 1]);
end

% make sure x direction is normal to make right/right
set(gca,'XDir','normal');
set(gca,'YDir','normal');
set(gca,'ZDir','normal');


%%%%%%%%%%%%%%%%%%%%%
%%   dispSurface   %%
%%%%%%%%%%%%%%%%%%%%%
function dispSurface

global gFlatViewer;
figure(gFlatViewer.f);

% get the patch vertices
patchVtcs = gFlatViewer.flat.patch2parent(:,2);

% get the vertexes/triangles and curvature
if gFlatViewer.whichSurface == 1
  vtcs = gFlatViewer.surfaces.outer.vtcs;
  tris = gFlatViewer.surfaces.outer.tris;
  c = gFlatViewer.curv;
elseif gFlatViewer.whichSurface == 2
  vtcs = gFlatViewer.surfaces.inner.vtcs;
  tris = gFlatViewer.surfaces.inner.tris;
  c = gFlatViewer.curv;
else
  vtcs = gFlatViewer.flat.vtcs;
  tris = gFlatViewer.flat.tris;
  c = gFlatViewer.curv(patchVtcs);
%  c = (c-min(c))./((max(c)-min(c)))>0.5;
  %  if this is a real flat, then view from above
  if ~isfield(gFlatViewer.flat,'radius')
    view([0 90]);
  end
end  

% clear the axis
cla;

% not sure why, but this is necessary to set up
% the axis so that right is right...
imagesc(0);

% get the colors that we want to show for that patch
[co alpha] = getPatchColoring;

% now set the overlay
if gFlatViewer.whichSurface <= 2
  overlay = NaN(length(c),3);
  overlay(patchVtcs,:) = co;
else
  overlay(:,:) = co;
end

% get the roi overlay
if isfield(gFlatViewer,'viewNum') && gFlatViewer.displayROIs
  % recompute roiOverlay
  if ~isfield(gFlatViewer,'roiOverlays') || ...
	length(gFlatViewer.roiOverlays) < gFlatViewer.whichSurface || ...
	isempty(gFlatViewer.roiOverlays{gFlatViewer.whichSurface})
    % get the vertices for which to calculate the roi overlay
    if gFlatViewer.whichSurface <= 2
      baseCoords = round(vtcs);
    else
      baseCoords = round(gFlatViewer.surfaces.inner.vtcs(patchVtcs,:));
    end
    % and compute them
    gFlatViewer.roiOverlays{gFlatViewer.whichSurface} = computeROIOverlay(baseCoords);
  end
  % get the overlay from the global
  roiOverlay = gFlatViewer.roiOverlays{gFlatViewer.whichSurface};
else
  roiOverlay = [];
end

% move vertices into center
vtcs(:,1) = vtcs(:,1)-mean(vtcs(:,1));
vtcs(:,2) = vtcs(:,2)-mean(vtcs(:,2));
vtcs(:,3) = vtcs(:,3)-mean(vtcs(:,3));

% convert the curvature to grayscale values
% note the 1.2 is because that is what we set
% the clim values to be to make them look nice
cmap = gray;
limval = 1.2;
c(c>limval) = limval;
c(c<-limval) = -limval;
c = round((size(cmap,1)-1)*(c-min(c))./(max(c)-min(c))+1);

% now make a combined overlay which has the grayscale
% values for the surface and the overlay values for
% where the patch is.
combinedOverlay(:,1) = cmap(c);
combinedOverlay(:,2) = cmap(c);
combinedOverlay(:,3) = cmap(c);
overlayPoints = ~isnan(overlay(:,1));
combinedOverlay(overlayPoints,:) = alpha*overlay(overlayPoints,:)+(1-alpha)*combinedOverlay(overlayPoints,:);

patch('vertices', vtcs, 'faces', tris,'FaceVertexCData', combinedOverlay, 'facecolor', 'interp','edgecolor', 'none');

% draw the surface and the overlay
%patch('vertices', vtcs, 'faces', tris, ...
%      'FaceVertexCData', c, ...
%      'facecolor', 'interp', ...
%      'edgecolor', 'none');
%patch('vertices', vtcs, 'faces', tris, ...
%      'FaceVertexCData', overlay, ...
%      'FaceColor', 'interp', 'Edgecolor','none','FaceAlpha',alpha);
if ~isempty(roiOverlay)
  patch('vertices', vtcs, 'faces', tris, ...
	'FaceVertexCData', roiOverlay, ...
	'FaceColor', 'interp', 'Edgecolor','none','FaceAlpha',.4);
end

% set axis stuff
axis off;axis equal;colormap(gray);axis tight;
camup('manual');
set(gca,'CLim',[-1.2 1.2]);

if gFlatViewer.whichSurface <= 2
  hPos = round(get(gFlatViewer.hSliders.h,'Value'));
  vPos = round(get(gFlatViewer.hSliders.v,'Value'));
  setViewAngle(hPos,vPos);
end
%%%%%%%%%%%%%%
% dispVolume
%%%%%%%%%%%%%%
function dispVolume(sliceIndex,slice)

global gFlatViewer;
figure(gFlatViewer.f);
cla(gca(gFlatViewer.f),'reset');

if length(size(gFlatViewer.anat.data)) < 3
  disp(sprintf('(mrFlatViewer) Could not display image becuase it does not have 3 dimensions'));
  return
end

% display a slice of the anatomy image
switch sliceIndex
  case {1}
   img = gFlatViewer.anat.data(slice,:);
  case {2}
   img = gFlatViewer.anat.data(:,slice,:);
  case {3}
   img = gFlatViewer.anat.data(:,:,slice);
end
imagesc(img);
colormap(gray);

axis image;
% axis off;
hold on

if min(img(:)) ~= max(img(:))
  set(gca,'CLim',[min(img(:)) max(img(:))]);
end
% display patch and white matter/gray matter
whichInx = gFlatViewer.flat.patch2parent(:,2);
wmPatchNodes = gFlatViewer.surfaces.inner.vtcs(whichInx,:);
gmPatchNodes = gFlatViewer.surfaces.outer.vtcs(whichInx,:);

% get full white matter/gray matter nodes
wmNodes = gFlatViewer.surfaces.inner.vtcs;
gmNodes = gFlatViewer.surfaces.outer.vtcs;

% Plot the nodes for the gray/white matter surfaces
wmNodes = wmNodes( find( round(wmNodes(:,sliceIndex))==slice), : );
plot(wmNodes(:,2), wmNodes(:,1), 'w.', 'markersize', 1);

gmNodes = gmNodes( find( round(gmNodes(:,sliceIndex))==slice), : );
plot(gmNodes(:,2), gmNodes(:,1), 'y.', 'markersize', 1);

% plot the patch nodes, displaying both deep and superficial surfaces
co = getPatchColoring;
if ~isnan(co(1))
  wmco = co(find( round(wmPatchNodes(:,sliceIndex))==slice),:);
  % make into magenta vs blue
  i = wmco(:,1);
  wmco(:,1) = i;
  wmco(:,2) = 0;
  wmco(:,3) = max(i,1-i);
else
  wmco(1:length(find( round(wmPatchNodes(:,sliceIndex))==slice)),1)=1;
  wmco(1:length(find( round(wmPatchNodes(:,sliceIndex))==slice)),2)=1;
  wmco(1:length(find( round(wmPatchNodes(:,sliceIndex))==slice)),3)=1;
end

wmPatchNodes = wmPatchNodes( find( round(wmPatchNodes(:,sliceIndex))==slice), : );

if any(gFlatViewer.patchColoring == [1 length(gFlatViewer.patchColoringTypes)])
  % draw all the points in the same color (if there are any)
  if ~ieNotDefined('wmco')
    plot(wmPatchNodes(:,2), wmPatchNodes(:,1), '.', 'markersize', 1,'Color',wmco(1,:));
  end
  % otherwise each pixel has to be set
else
  for i = 1:length(wmPatchNodes(:,1))
    plot(wmPatchNodes(i,2), wmPatchNodes(i,1), '.', 'markersize', 1,'Color',wmco(i,:)');
  end
end

if ~isnan(co(1))
  gmco = co(find( round(gmPatchNodes(:,sliceIndex))==slice),:);
else
  gmco(1:length(find( round(gmPatchNodes(:,sliceIndex))==slice)),1)=1;
  gmco(1:length(find( round(gmPatchNodes(:,sliceIndex))==slice)),2)=1;
  gmco(1:length(find( round(gmPatchNodes(:,sliceIndex))==slice)),3)=0;
end
gmPatchNodes = gmPatchNodes( find( round(gmPatchNodes(:,sliceIndex))==slice), : );
% uniform patch coloring
if any(gFlatViewer.patchColoring == [1 length(gFlatViewer.patchColoringTypes)])
  % draw all the points in the same color (if there are any)
  if ~ieNotDefined('gmco')
    plot(gmPatchNodes(:,2), gmPatchNodes(:,1), '.', 'markersize', 1,'Color',gmco(1,:));
  end
  % otherwise each pixel has to be set
else
  for i = 1:length(gmPatchNodes(:,1))
    plot(gmPatchNodes(i,2), gmPatchNodes(i,1), '.', 'markersize', 1,'Color',gmco(i,:));
  end
end

view([0 90]);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   getPatchColoring   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [co alpha] = getPatchColoring

global gFlatViewer;

alpha = 0.6;

% get the patch vertices
patchVtcs = gFlatViewer.flat.patch2parent(:,2);

% one is uniform, 2-4 are red/blue
switch gFlatViewer.patchColoring
 % uniform
 case 1
  % make everybody red
  co = ones(1,gFlatViewer.flat.Nvtcs);
  % anatomical directions
 case {2,3,4}
  co = gFlatViewer.surfaces.outer.vtcs(patchVtcs,gFlatViewer.patchColoring-1)';
  co = (co-min(co))./(max(co)-min(co));
  % curvature
 case {5,6}
  % get curvature
  curv = gFlatViewer.curv(patchVtcs)';
  curv = (curv-min(curv))./(max(curv)-min(curv));
  if gFlatViewer.patchColoring == 6
    curv = 1-curv;
  end
  % flatten distribution
  co = flattenDistribution(curv);
 % Areal distortion
 case {7,8,9,10}
  tris = gFlatViewer.flat.tris;
  % get area in volume
  if gFlatViewer.patchColoring == 9
    % use outer surface
    volumeArea = getTriangleArea(gFlatViewer.surfaces.outer.vtcs(patchVtcs,:),tris);
  else
    % use inner surface
    volumeArea = getTriangleArea(gFlatViewer.surfaces.inner.vtcs(patchVtcs,:),tris);
  end
  patchArea = getTriangleArea(gFlatViewer.flat.vtcs,tris);
  % get the distortion as the ratio of the area in the
  % patch to the volume
  trisDistortion = patchArea./volumeArea;
  % now convert this into a color for each vertex
  distortion = ones(1,length(patchVtcs));
  % note that this is a shortcut, it just
  % sets each vertex to one value of the distortion
  % even though the vertex may belong to many triangles
  distortion(tris(:,1)) = trisDistortion;
  distortion(tris(:,2)) = trisDistortion;
  distortion(tris(:,3)) = trisDistortion;
  % normalize
  if gFlatViewer.patchColoring < 9
    co = (distortion-min(distortion))./(max(distortion)-min(distortion));
    co = flattenDistribution(co);
    % invert colors
    distortion = log10(distortion);
    if gFlatViewer.patchColoring == 8
      co = 1-co;
      title(sprintf('Stretch: max=%0.2fx median=%0.2fx mean=%0.2fx',10^max(distortion),10^median(distortion(distortion>0)),10^mean(distortion(distortion>0))));
    else
      title(sprintf('Compression: max=%0.2fx median=%0.2fx mean=%0.2fx',10^abs(min(distortion)),10^abs(median(distortion(distortion<0))),10^abs(mean(distortion(distortion<0)))));
    end
  else
    % set the colors to the absolute value
    % of the log of the distortion. This scales
    % so that doubling or halving the area from
    % the volume to the patch give you the same
    % number. Anything above 10x distortion is 1
    distortion = abs(log10(distortion));
    co = distortion;
    co(co>1) = 1;
    % print out some statistics
    distortion = 10.^distortion;
    title(sprintf('Distortion: max=%0.2fx median=%0.2fx mean=%0.2fx',max(distortion),median(distortion),mean(distortion)));
  end
 %current overlay
 case {11,12}
  whichInx = gFlatViewer.flat.patch2parent(:,2);
  % get the coordinates from the right surface
  if gFlatViewer.whichSurface == 1
    baseCoords = gFlatViewer.surfaces.outer.vtcs(whichInx,:);
  else
    baseCoords = gFlatViewer.surfaces.inner.vtcs(whichInx,:);    
  end
  % make homogenous
  baseCoords(:,4) = 1;
  baseCoords = baseCoords';
  baseDims = [size(baseCoords,2) 1];
  % get the view
  v = viewGet([],'view',gFlatViewer.viewNum);
  if ~isempty(v) & ~isempty(viewGet(v,'currentOverlay'))
    % get the base2scan xform
    base2scan = viewGet(v,'base2scan');
    % and get the overlay
    overlay = computeOverlay(v,base2scan,baseCoords,baseDims);
    overlay.RGB(overlay.alphaMap==0) = nan;
    co = squeeze(overlay.RGB);
    alpha = viewGet(v,'alpha');
  else
    co = zeros(size(baseCoords,2),3);
    co(:) = nan;
    alpha = 1;
  end
  if gFlatViewer.patchColoring==12
    % get curvature
    curv = gFlatViewer.curv(patchVtcs)';
    % the default curvature mapping is to use the
    % gray colortable clipping values at -1.2 and 1.2
    % and mapping them lineraly into the colortable.
    % This is what is specified by sending a scalar value
    % to the patch command. So we replicate that here to
    % draw the patch
    cmap = brighten(gray,-0.5);
    curv(curv>1.2) = 1.2;
    curv(curv<-1.2) = -1.2;
    curv = round((size(cmap,1)-1)*(curv-min(curv))./(max(curv)-min(curv))+1);
    % set all the points that don't show up in the overlay
    % with the curvature slightly darkened
    noOverlayPoints = find(isnan(co(:,1)));
    co(noOverlayPoints,1) = cmap(curv(noOverlayPoints),1);
    co(noOverlayPoints,2) = cmap(curv(noOverlayPoints),2);
    co(noOverlayPoints,3) = cmap(curv(noOverlayPoints),3);
    % set all the points that are in the overlay to be correctply
    % alpha blended with the patch points
    overlayPoints = find(~isnan(co(:,1)));
    co(overlayPoints,1) = alpha*co(overlayPoints,1)+(1-alpha)*cmap(curv(overlayPoints),1);
    co(overlayPoints,2) = alpha*co(overlayPoints,2)+(1-alpha)*cmap(curv(overlayPoints),2);
    co(overlayPoints,3) = alpha*co(overlayPoints,3)+(1-alpha)*cmap(curv(overlayPoints),3);
    % now set the alpha to 1, since we have already done the alpha
    % blending here.
    alpha = 1;
  end
  return
 % no coloring
 case {length(gFlatViewer.patchColoringTypes)}
  % make everyone nan
  co = ones(gFlatViewer.flat.Nvtcs,3);
  co(:) = nan;
  return
end

% intermediate values turn to gray, so avoid them
co((co>0.3)&(co<0.5)) = 0.3;
co((co>0.5)&(co<0.7)) = 0.7;

% make into RGB
co(2:3,:) = [1-co;1-co];
co = co';

%%%%%%%%%%%%%%%%%%%%%%%%%
%%   getTriangleArea   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function area = getTriangleArea(vtcs,tris)

% get length of a each side of triangles
a = sqrt(sum((vtcs(tris(:,1),:)-vtcs(tris(:,2),:))'.^2));
b = sqrt(sum((vtcs(tris(:,2),:)-vtcs(tris(:,3),:))'.^2));
c = sqrt(sum((vtcs(tris(:,3),:)-vtcs(tris(:,1),:))'.^2));

% get semiperimeter (i.e. 1/2 perimeter)
p = (a+b+c)/2;

% use Heron's formula for size of triangle given side lengths
area = sqrt(p.*(p-a).*(p-b).*(p-c));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   flattenDistribution   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function co = flattenDistribution(co,nbins)

if ~exist('nbins','var'),nbins = 10;end

% sort the values
[co sortIndex] = sort(co);
% get binsize
binSize = floor(length(co)/nbins);
% add even number of values back into each bin
for i = 1:nbins
  co(sortIndex((i-1)*binSize+1:min(length(co),i*binSize))) = i/nbins;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   computeROIOverlay   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function roiOverlay = computeROIOverlay(baseCoords);

global gFlatViewer;
roiOverlay = [];
disppercent(-inf,'(mrFlatViewer) Computing ROI Overlay');

% get view information
v = viewGet([],'view',gFlatViewer.viewNum);
numROIs = viewGet(v,'numROIs');
baseVoxelSize = [1 1 1];
  
% init the overlay
roiOverlay = zeros(size(baseCoords,1),3);
roiOverlay(:) = nan;

% deal with selected ROI color
selectedROI = viewGet(v,'currentroi');
selectedROIColor = mrGetPref('selectedROIColor');

% get which ROIs to do
showROIs = viewGet(v,'showROIs');
if strcmp(showROIs,'none')
  return
end
if strfind(showROIs,'selected')
  rois = selectedROI;
else
  rois = 1:numROIs;
end

for roinum = rois
  % get ROI info
  roiCoords = viewGet(v,'roiCoords',roinum);
  base2roi = viewGet(v,'base2roi',roinum);
  roiVoxelSize = viewGet(v,'roiVoxelSize',roinum);
  if roinum ~= selectedROI
    roiColorRGB = viewGet(v,'roiColorRGB',roinum);
  else
    if strcmp(selectedROIColor,'none')
      roiColorRGB = viewGet(v,'roiColorRGB',roinum);
    else
      roiColorRGB = color2RGB(selectedROIColor);
    end
  end
  % get the base coord that match the roi
  roiBaseCoords = round(xformROIcoords(roiCoords,inv(base2roi),roiVoxelSize,baseVoxelSize));
  roiBaseCoords = roiBaseCoords(1:3,:)';
  roiVertices = find(ismember(baseCoords,roiBaseCoords,'rows'));
  % and set them to the roi color
  roiOverlay(roiVertices,1) = roiColorRGB(1);
  roiOverlay(roiVertices,2) = roiColorRGB(2);
  roiOverlay(roiVertices,3) = roiColorRGB(3);
  disppercent(roinum/numROIs);
end
disppercent(inf);

%%%%%%%%%%%%%%%%%%%%%%%
%%   switchAnatomy   %%
%%%%%%%%%%%%%%%%%%%%%%%
function switchAnatomy(params)

global gFlatViewer;

% load the anatomy and view
disppercent(-inf,sprintf('(mrFlatViewer) Load %s',params.anatFileName));
[gFlatViewer.anat.data gFlatViewer.anat.hdr] = mlrImageReadNifti(fullfile(params.path, params.anatFileName));
gFlatViewer = xformSurfaces(gFlatViewer);
% switch to 3D anatomy view
global gParams
gFlatViewer.whichSurface = 3;
set(gParams.ui.varentry{1},'Value',gFlatViewer.whichSurface)
refreshFlatViewer([],[],1);
disppercent(inf);

%%%%%%%%%%%%%%%%%%%%
%%   switchFlat   %%
%%%%%%%%%%%%%%%%%%%%
function switchFlat(params)

global gFlatViewer;

% load the anatomy and view
disppercent(-inf,sprintf('(mrFlatViewer) Load %s',params.flatFileName));
gFlatViewer.flat = loadSurfOFF(fullfile(params.path, params.flatFileName));
% switch to flat view
global gParams
refreshFlatViewer([],[],1);
disppercent(inf);

%%%%%%%%%%%%%%%%%%%%
%%   switchFile   %%
%%%%%%%%%%%%%%%%%%%%
function switchFile(whichSurface,params)

global gFlatViewer;

% if the user wants to find a new file
addFilename = 0;
if strcmp(params.(whichSurface),'Find file')
  if strcmp(whichSurface,'curvFileName')
    [filename, pathname] = uigetfile({'*.vff','VFF Curvature files (*.vff)'},'Select curvature file',gFlatViewer.path);
    whichControl = gFlatViewer.guiloc.filenames+3;
  else
    [filename, pathname] = uigetfile({'*.off','OFF Surface files (*.off)'},'Select surface',gFlatViewer.path);
    whichControl = gFlatViewer.guiloc.filenames+1+find(strcmp(whichSurface,{'outerCoordsFileName','innerCoordsFileName'}));
  end
  filename = getRelativePath(gFlatViewer.path,fullfile(pathname,filename));
  addFilename = 1;
else
  filename = params.(whichSurface);
end

% try to load it
disppercent(-inf,sprintf('(mrFlatViewer) Loading %s',filename));
if filename ~= 0
  if strcmp(whichSurface,'curvFileName')
    file = myLoadCurvature(fullfile(params.path, filename));
    whichControl = gFlatViewer.guiloc.filenames+3;;
  else
    file = myLoadSurface(fullfile(params.path, filename));
    whichControl = gFlatViewer.guiloc.filenames+1+find(strcmp(whichSurface,{'outerCoordsFileName','innerCoordsFileName'}));
  end
else
  file = [];
end

% get the proper field name.
whichSurfaceTypes = {'outerCoordsFileName','innerCoordsFileName','curvFileName','anatFileName'};
whichFieldName    = {'outer',              'inner',              'curv',        'anat'        };
surfaceFieldName = whichFieldName{find(strcmp(whichSurface,whichSurfaceTypes))};

% get which surface field it is
if ~isempty(file)
  if strcmp(whichSurface,'curvFileName')
    gFlatViewer.(surfaceFieldName)=file;
  else
    gFlatViewer.surfaces.(surfaceFieldName)=file;
    gFlatViewer = xformSurfaces(gFlatViewer);
    % set the correct one to display
    gFlatViewer.whichSurface = find(strcmp(whichSurface,{'outerCoordsFileName','innerCoordsFileName'}));
  end    
  % and change the ui control
  global gParams;
  set(gParams.ui.varentry{gFlatViewer.guiloc.whichSurface},'Value',gFlatViewer.whichSurface)
  % add the filename to the control if necessary
  if addFilename
    currentChoices = get(gParams.ui.varentry{gFlatViewer.whichSurface+1},'String');
    currentChoices = setdiff(currentChoices,'Find file');
    currentChoices = putOnTopOfList(filename,currentChoices);
    currentChoices{end+1} = 'Find file';
    set(gParams.ui.varentry{gFlatViewer.guiloc.filenames+gFlatViewer.whichSurface},'String',currentChoices)
    set(gParams.ui.varentry{gFlatViewer.guiloc.filenames+gFlatViewer.whichSurface},'Value',1)
  end
else
  global gParams;
  % switch back to first on list
  currentChoices = get(gParams.ui.varentry{whichControl},'String');
  set(gParams.ui.varentry{whichControl},'Value',1)
  if ~strcmp(whichSurface,'curv')
    gFlatViewer.surfaces.(whichSurface) = myLoadSurface(fullfile(params.path, currentChoices{1}));
  else
    gFlatViewer.curv = myLoadCurvature(fullfile(params.path, currentChoices{1}));
  end
end
refreshFlatViewer([],[],1);

disppercent(inf);

%%%%%%%%%%%%%%%%%%%%%%%
%%   myLoadSurface   %%
%%%%%%%%%%%%%%%%%%%%%%%
function surf = myLoadSurface(filename,filepath,onlyLoadHeader)

% default to loading data
if ieNotDefined('onlyLoadHeader'),onlyLoadHeader = 0;end

% put on path
if ~ieNotDefined('filepath')
  filename = fullfile(filepath,filename);
end

global gFlatViewer;
% load the surface
surf = loadSurfOFF(filename,onlyLoadHeader);
if isempty(surf),return,end

% check that it has the correct number of vertices
%if ~isequal(gFlatViewer.flat.nParent(1:2),[surf.Nvtcs surf.Ntris]')
if ~isequal(gFlatViewer.flat.nParent(1),surf.Nvtcs);
  % dispaly warning, but only if mismatchWarning is set,
  % this way when we first load surfaces just for checking it
  % won't complain
  if gFlatViewer.mismatchWarning
    mrWarnDlg(sprintf('(mrFlatViewer) Surface %s does not match patch Nvtcs: %i vs %i, Ntris: %i vs %i',filename,gFlatViewer.flat.nParent(1),surf.Nvtcs,gFlatViewer.flat.nParent(2),surf.Ntris));
  end
  surf = [];
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%   myLoadCurvature   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function curv = myLoadCurvature(filename,filepath,onlyLoadHeader)

% default to loading data
if ieNotDefined('onlyLoadHeader'),onlyLoadHeader = 0;end

% put on path
if ~ieNotDefined('filepath')
  filename = fullfile(filepath,filename);
end

global gFlatViewer;
% load the curvature
[curv curvhdr] = loadVFF(filename,onlyLoadHeader);
if isempty(curvhdr)
  return
end

% % check that it has the correct number of vertices
if ~isequal(gFlatViewer.flat.nParent(1), curvhdr.size(3))
  % dispaly warning, but only if mismatchWarning is set,
  % this way when we first load surfaces just for checking it
  % won't complain
  if gFlatViewer.mismatchWarning
    mrWarnDlg(sprintf('(mrFlatViewer) Curvature file %s does not match patch Nvtcs: %i vs %i',filename,gFlatViewer.flat.nParent(1), curvhdr.size(3)));
  end
  curv = [];
  return
end

% return header only
if onlyLoadHeader
  curv = curvhdr;
else 
  curv = curv';
end

% function that handles conversion of surface vtcs from 
% Jonas' world coordinates to array coordinates. This
% is usually just an offset that is read from the nifti
% header as the undocumented fields pixdim(6:8), but
% may in the future be a complete nifti style xform
function gFlatViewer = xformSurfaces(gFlatViewer)

surfaces = {'inner','outer'};
for surfNum = 1:length(surfaces)
  % and store them back as the vtcs
  gFlatViewer.surfaces.(surfaces{surfNum}) = xformSurfaceWorld2Array(gFlatViewer.surfaces.(surfaces{surfNum}),gFlatViewer.anat.hdr);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   makeFlatFromRadius   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function flat = makeFlatFromRadius(flat,radius,startPoint,surfFileName)

if nargin == 4
  % make sure we have an off
  surfFileName = setext(surfFileName,'off');
  % get the surface name and load it
  flat.parentSurfaceName = getLastDir(surfFileName);
  surf = loadSurfOFF(fullfile(flat.path,surfFileName));
  if isempty(surf),flat =[];return;end
  % create a connection matrix
  mesh.uniqueVertices = surf.vtcs;
  mesh.uniqueFaceIndexList = surf.tris;
  mesh.connectionMatrix = findConnectionMatrix(mesh);
  % note that here, we could pass in a scaling. As long
  % as the volume is 1x1x1 mm, the scaling is in mm though.
  flat.distanceMatrix = find3DNeighbourDists(mesh);
  % keep the parent name
  flat.parentSurfaceName = surfFileName;
  % and parent info
  flat.nParent = [surf.Nvtcs surf.Ntris surf.Nedges];
  flat.parent.tris = surf.tris;
  flat.parent.vtcs = surf.vtcs;
end
if nargin >= 3
  % get the nearest vertex to the start point in the surface
  flat.startVertex = assignToNearest(flat.parent.vtcs,startPoint);
  % and remember those coordinates
  flat.startPoint = round(flat.parent.vtcs(flat.startVertex,:));
  % compute distance from start vertex to every other vertex
  flat.distance = dijkstra(flat.distanceMatrix,flat.startVertex);
end
flat.radius = radius;
% now get what vertexes are within specified radius
flat.patch2parent = [];
flat.patch2parent(:,1) = 1:sum(flat.distance < flat.radius);
flat.patch2parent(:,2) = find(flat.distance < flat.radius);
% fill in some fields
flat.Nvtcs = size(flat.patch2parent,1);
% get the triangles corresponding to this patch
whichParentTris = ismember(flat.parent.tris,flat.patch2parent(:,2));
whichParentTris = find(sum(whichParentTris')==3);
flat.tris = flat.parent.tris(whichParentTris,:);
flat.Ntris = size(flat.tris,1);
% convert those tris to flat vertexs
[tf flat.tris] = ismember(flat.tris(:),flat.patch2parent(:,2));
flat.tris = reshape(flat.tris,flat.Ntris,3);
% get vertices
flat.vtcs = flat.parent.vtcs(flat.patch2parent(:,2),:);
% fill out rest of fields
flat.Nedges = flat.Nvtcs+flat.Ntris-1;
flat.nPatch = [flat.Nvtcs flat.Ntris flat.Nedges];
% set the name of the patch
flat.name = sprintf('%s_Flat_%i_%i_%i_Rad%i.off',stripext(flat.parentSurfaceName),flat.startPoint(1),flat.startPoint(2),flat.startPoint(3),flat.radius);

%%%%%%%%%%%%%%%%%%%%%%%
%%   setFlatRadius   %%
%%%%%%%%%%%%%%%%%%%%%%%
function setFlatRadius(params)

global gFlatViewer;

% see if we need to change name
updateName = 0;
if isempty(params.flatFileName) || strcmp(gFlatViewer.flat.name,params.flatFileName)
  updateName = 1;
end

% reset the patch for the current selected radius
gFlatViewer.flat = makeFlatFromRadius(gFlatViewer.flat,params.radius);

% and update name
if (updateName)
  params.flatFileName = gFlatViewer.flat.name;
end

% reset parameters
mrParamsSet(params);

% and refresh
refreshFlatViewer([],[],1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   setFlatStartPoint   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setFlatStartPoint(params)

global gFlatViewer;

% see if we need to change name
updateName = 0;
if isempty(params.flatFileName) || strcmp(gFlatViewer.flat.name,params.flatFileName)
  updateName = 1;
end

% reset the patch for the current selected radius
gFlatViewer.flat = makeFlatFromRadius(gFlatViewer.flat,params.radius,[params.x params.y params.z]);

% and update location of startPoint
params.x = gFlatViewer.flat.startPoint(1);
params.y = gFlatViewer.flat.startPoint(2);
params.z = gFlatViewer.flat.startPoint(3);

% and update name
if (updateName)
  params.flatFileName = gFlatViewer.flat.name;
end

% reset parameters
mrParamsSet(params);

% and refresh
refreshFlatViewer([],[],1);
