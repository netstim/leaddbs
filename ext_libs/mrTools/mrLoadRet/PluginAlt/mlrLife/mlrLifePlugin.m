% mlrLifePlugin
%
%        $Id:$ 
%      usage: mlrLifePlugin(action,<v>)
%         by: justin gardner & franco pestilli
%       date: 09/09/2014
%    purpose: Plugin function for LiFE
%
function retval = mlrLifePlugin(action,v)

% check arguments
if ~any(nargin == [1 2])
  help DefaultPlugin
  return
end

switch action
 case {'install','i'}
  % check for a valid view
  if (nargin ~= 2) || ~isview(v)
     disp(sprintf('(mlrLifePlugin) Need a valid view to install plugin'));
  else
    % add the menu item for Diffusion
    mlrAdjustGUI(v,'add','menu','Diffusion','/File/ROI');
    % add one to open Dt6 file
    mlrAdjustGUI(v,'add','menu','Load Dt6','/File/Diffusion/','Callback',@mlrLifeLoadDt6);
    % this menu item allows importation of fascicles as a surface
    mlrAdjustGUI(v,'add','menu','Import fascicles','/File/Base anatomy/Import surface','Callback',@mlrLifeImportFascicles);

    % return true to indicate successful plugin
    retval = true;
   end
 % return a help string
 case {'help','h','?'}
   retval = 'This is an example plugin, it just installs a menu item to Select Plugins.';
 otherwise
   disp(sprintf('(mlrLifePlugin) Unknown command %s',action));
end

%%%%%%%%%%%%%%%%%%%%%%%
%    mlrLifeLoadDt6   %
%%%%%%%%%%%%%%%%%%%%%%%
function mlrLifeLoadDt6(hObject,eventdata)

% code-snippet to get the view from the hObject variable. Not needed for this callback.
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% check system
if ~mlrLifeSystemCheck,return,end

% bring up file dialog to select dt6 file
global mlrLifePluginStartPathStr;
if isempty(mlrLifePluginStartPathStr)
  mlrLifePluginStartPathStr = pwd;
end
filterspec = {'dt6.mat','Diffusion tensor 6 parameter file';'*.*','All files'};
title = 'Choose diffusion tensor 6 parameter file';
dt6Filename = mlrGetPathStrDialog(mlrLifePluginStartPathStr,title,filterspec,'off');
if isempty(dt6Filename),return,end

% use mrDiffusion functions to load dti
[dt, t1, o] = dtiLoadDt6(dt6Filename);
if isempty(dt),return,end

% save the path
mlrLifePluginStartPathStr = fileparts(dt6Filename);

% now bring up selection dialog
paramsInfo = {{'groupName','Diffusion'},...
	      {'analysisName','Dt6'}};
params = mrParamsDialog(paramsInfo);
if isempty(params),return,end

% convert t1 into a MLR anatomy
% Set required structure fields (additional fields are set to default
% values when viewSet calls isbase).
base.name = stripext(stripext(getLastDir(dt.files.t1)));
base.data = double(t1.img);
% create a header
h.dim = size(t1.img);
[tf h] = mlrImageIsHeader(h);
h.pixdim = t1.mmPerVoxel;
h.qform = diag([h.pixdim 1]);
h.sform = t1.xformToAcpc;

% set nifti header and permutation matrix in base structure
base.hdr = mlrImageGetNiftiHeader(h);
base.permutationMatrix = getPermutationMatrix(base.hdr);

% make it a full base
[tf base] = isbase(base);

% add it to the view
v = viewSet(v,'newBase',base);

% make group
v = viewSet(v,'newGroup',params.groupName);
v = viewSet(v,'curGroup',params.groupName);

% make nifti header for DTI
h = [];hdr = [];
h.dim = size(dt.dt6);
[tf h] = mlrImageIsHeader(h);
h.pixdim = dt.mmPerVoxel;
h.qform = diag([h.pixdim 1]);
h.sform = dt.xformToAcpc;
hdr = mlrImageGetNiftiHeader(h);

% save tseries
saveNewTSeries(v,dt.dt6,[],hdr);
scanNum = viewGet(v,'nScans');
v = viewSet(v,'curScan',scanNum);

% make the analysis
a = mlrMakeAnalysis(v,params.analysisName);
% add the b0 Overlay
a = mlrMakeAnalysis(v,a,double(o.b0),'overlayName','b0','overlayColormap',hot(312));

% convert RGB image into an indexed image
cmap = [gray(16);hsv(48);hsv(48)*0.8;hsv(48)*0.6;hsv(48)*0.4;hsv(48)*0.2];
RGB = reshape(o.vectorRGB,prod(h.dim(1:2)),h.dim(3),3);
[RGB cmap] = rgb2ind(RGB,cmap,'nodither');
RGB = reshape(RGB,h.dim(1:3));
% add it as an overlay
a = mlrMakeAnalysis(v,a,double(RGB)/size(cmap,1),'overlayName','vectorRGB','overlayColormap',cmap,'overlayRange',[0 1]);

% set the analysis in the view
v = viewSet(v,'newAnalysis',a);
% save it
saveAnalysis(v,a.name);
% and refresh the display
refreshMLRDisplay(v);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrLifeImportFascicles    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrLifeImportFascicles(hObject,eventdata)

% code-snippet to get the view from the hObject variable. Not needed for this callback.
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% check system
if ~mlrLifeSystemCheck,return,end

% Load the fascicles from disk
[fg, fgFilename] = fgRead;
if isempty(fg), return,end

% bring up file dialog to select dt6 file
% this is so that we can get the T1 file
% FIX, FIX, FIX - should find some way to automatically
% load the T1 file without having the user have to go find it.
startPathStr = fullfile(fileparts(fgFilename),'..');
filterspec = {'dt6.mat','Diffusion tensor 6 parameter file';'*.*','All files'};
title = 'Choose diffusion tensor 6 parameter file';
dt6Filename = mlrGetPathStrDialog(startPathStr,title,filterspec,'off');

% use mrDiffusion functions to load dti
[dt, t1, o] = dtiLoadDt6(dt6Filename);
if isempty(dt),return,end

% Get the header from the T1 so that we have the correct sform etc
h.dim = size(t1.img);
[tf h] = mlrImageIsHeader(h);
h.pixdim = t1.mmPerVoxel;
h.qform = diag([h.pixdim 1]);
h.sform = t1.xformToAcpc;

% convert the X, Y, Z coordinates which are in AC/PC back to 3D image coordinates
huh = h.sform;
h.sform = eye(4);
%huh(1,1) = 1*t1.mmPerVoxel(1);
%huh(2,2) = 1*t1.mmPerVoxel(2);
%huh(3,3) = 1*t1.mmPerVoxel(3);
%huh(1,4) = huh(1,4)*t1.mmPerVoxel(1);
%huh(2,4) = huh(2,4)*t1.mmPerVoxel(2);
%huh(3,4) = huh(3,4)*t1.mmPerVoxel(3);
%scaleXform = diag([t1.mmPerVoxel 1]);
%fg = dtiXformFiberCoords(fg,scaleXform);
%huh = huh*inv(scaleXform);
%fg = dtiXformFiberCoords(fg,inv(huh));

% Build frames from the fascicles
[X, Y, Z] = mbaBuildFascicleFrame(fg.fibers);
%[Y, X, Z] = mbaBuildFascicleFrame(fg.fibers);

% number of fascicles
nFascicles = length(X);

% Build a patch from the frame
nTotalVertices = 0;nTotalTris = 0;
for i = 1:nFascicles
  fasciclePatches{i} = surf2patch(X{i},Y{i},Z{i},'triangles');
  % compute how many vertices and tris we have all together
  nTotalVertices = size(fasciclePatches{i}.vertices,1) + nTotalVertices;
  nTotalTris = size(fasciclePatches{i}.faces,1) + nTotalTris;
end

% create an MLR base structure
fascicleBase.name = fixBadChars(fg.name);

% set this to be a surface
fascicleBase.type = 2; 

% set nifti header and permutation matrix in base structure
fascicleBase.hdr = mlrImageGetNiftiHeader(h);
fascicleBase.permutationMatrix = getPermutationMatrix(fascicleBase.hdr);
fascicleBase.vol2mag = fascicleBase.hdr.sform44;

% set path and names of files
fascicleBase.coordMap.path = fileparts(fgFilename);
fascicleBase.coordMap.innerSurfaceFileName = getLastDir(fgFilename);
fascicleBase.coordMap.outerSurfaceFileName = getLastDir(fgFilename);
fascicleBase.coordMap.innerCoordsFileName = getLastDir(fgFilename);
fascicleBase.coordMap.outerCoordsFileName = getLastDir(fgFilename);
fascicleBase.coordMap.curvFileName = '';
fascicleBase.coordMap.anatFileName = dt.files.t1;

% initalize fields
fascicleBase.data = [];
% set dimensions for coordMap
fascicleBase.coordMap.dims = h.dim;
fascicleBase.coordMap.innerCoords = nan(1,nTotalVertices,1,3);
fascicleBase.coordMap.innerVtcs = nan(nTotalVertices,3);
fascicleBase.coordMap.tris = nan(nTotalTris,3);
nRunningTotalVertices = 0;
nRunningTotalTris = 0;

% now put all fascicles vertices and triangles into one coordMap
disppercent(-inf,sprintf('(mlrLifePlugin) Converting %i fascicles',nFascicles));
for iFascicle = 1:nFascicles
  % number of vertices and triangles
  nVertices = size(fasciclePatches{iFascicle}.vertices,1);
  nTris = size(fasciclePatches{iFascicle}.faces,1);
  % the data which is the grayscale value to color the fascicles with (rand for now)
  fascicleBase.data = [fascicleBase.data rand(1,nVertices)];
  % convert vertices to a coord map which has one x,y,z element for each possible
  % location on the surface (which actually is just a 1xnVerticesx1 image)
  % add these vertices to existing vertices
  fascicleBase.coordMap.innerCoords(1,nRunningTotalVertices+1:nRunningTotalVertices+nVertices,1,1) = fasciclePatches{iFascicle}.vertices(:,1);
  fascicleBase.coordMap.innerCoords(1,nRunningTotalVertices+1:nRunningTotalVertices+nVertices,1,2) = fasciclePatches{iFascicle}.vertices(:,2);
  fascicleBase.coordMap.innerCoords(1,nRunningTotalVertices+1:nRunningTotalVertices+nVertices,1,3) = fasciclePatches{iFascicle}.vertices(:,3);
  % these are the display vertices which are the same as the coords
  fascicleBase.coordMap.innerVtcs(nRunningTotalVertices+1:nRunningTotalVertices+nVertices,:) = fasciclePatches{iFascicle}.vertices;
  % triangle faces
  fascicleBase.coordMap.tris(nRunningTotalTris+1:nRunningTotalTris+nTris,:) = (fasciclePatches{iFascicle}.faces + nRunningTotalVertices);
  % update runing totals
  nRunningTotalVertices = nRunningTotalVertices + nVertices;
  nRunningTotalTris= nRunningTotalTris + nTris;
  disppercent(iFascicle/nFascicles);
end
disppercent(inf);

% copy the inner to outer since they are all the same for fascicles
fascicleBase.coordMap.outerCoords = fascicleBase.coordMap.innerCoords;
fascicleBase.coordMap.outerVtcs = fascicleBase.coordMap.innerVtcs;

% save the individual fascicles so we can rebuild later
fascicleBase.fascicles.patches = fasciclePatches;
fascicleBase.fascicles.n = length(fasciclePatches);
fascicleBase.fascicles.nTotalVertices = nTotalVertices;
fascicleBase.fascicles.nTotalTris = nTotalTris;

% make it a full base
[tf fascicleBase] = isbase(fascicleBase);

% add it to the view
v = viewSet(v,'newBase',fascicleBase);

% refresh the display to draw it
refreshMLRDisplay(v);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   mlrLifeSystemCheck   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = mlrLifeSystemCheck 

tf = false;
% check if we have fgRead from vistasoft
if exist('fgRead') ~= 2
  mrWarnDlg(sprintf('(mlrLifePlugin) You need to have vista installed to access function fgRead.\ngit clone https://github.com/vistalab/vistasoft.git vistasoft'));
  return
end

% check if we have mbaBuildFascicleFrame from mba
if exist('mbaBuildFascicleFrame') ~= 2
  mrWarnDlg(sprintf('(mlrLifePlugin) You need to have mba installed to access function mbaBuildFascicleFrame.\ngit clone https://github.com/francopestilli/mba.git mba'));
  return
end
 tf = true;