% importCaret.m
%
%        $Id$ 
%      usage: mlrImportCaret()
%         by: justin gardner
%       date: 12/03/09
%    purpose: function to import caret surfaces should be used after running mypreborder/postborder.
%             Should be run from within a FreeSurfer directory. mlrImportFreeSurfer should have
%             already been run on that directory.
%    options: verbose=0 :set to 1 to for more comments
%             atlasDir='../PALS_B12.LR': The directory where the atlas and the deformed surfaces
%                        for the subject should be found
%             doCurvature=1: Set to 0 to not compute curvature at all. Set to 2 to force recomputing%                        of curvature
%
function retval = mlrImportCaret(varargin)

% check arguments
if ~any(nargin == [0 1])
  help importCaret
  return
end

verbose = [];atlasDir = [];doCurvature = [];
getArgs(varargin,{'verbose=0','atlasDir=../PALS_B12.LR','doCurvature=1'});

% we need files from these places
surfRelaxDir = 'surfRelax';
caretFileDir = 'surf';

% check to make sure the directories exist
if ~checkDirs(caretFileDir,surfRelaxDir,atlasDir),return,end

% look for coord files in caret directory. These will be used for computing
% the xform that goes from 711-2B back to the original coordinates
[rightCoordFiles d.rightStemName] = getCoordFiles(caretFileDir,'R.Midthickness');
if isempty(rightCoordFiles),return,end
[leftCoordFiles d.leftStemName] = getCoordFiles(caretFileDir,'L.Midthickness');
if isempty(leftCoordFiles),return,end

% now compute the transformation from caret 711-2B back to original coordinates;
d.leftXform = computeTransformBetweenSurfaces(leftCoordFiles{2},leftCoordFiles{1});
if isempty(d.leftXform),return,end
d.rightXform = computeTransformBetweenSurfaces(rightCoordFiles{2},rightCoordFiles{1});
if isempty(d.rightXform),return,end

% now look for the files in the atlas directory
[topo d.leftAtlasAligned] = listCaretDir(sprintf('dirname=%s',atlasDir),sprintf('matchStr=%s',d.leftStemName),'noDisplay=1');
if isempty(d.leftAtlasAligned),disp(sprintf('(mlrImportCaret) Could not find deformed surface for %s in %s',d.leftStemName,atlasDir));return,end
[topo d.rightAtlasAligned] = listCaretDir(sprintf('dirname=%s',atlasDir),sprintf('matchStr=%s',d.rightStemName),'noDisplay=1');
if isempty(d.rightAtlasAligned),disp(sprintf('(mlrImportCaret) Could not find deformed surface for %s in %s',d.rightStemName,atlasDir));return,end

% look for topo files for LEFT and RIGHT hemispheres
[d.leftAtlasTopo d.rightAtlasTopo] = getAtlasTopo(atlasDir);
if isempty(d.leftAtlasTopo) || isempty(d.rightAtlasTopo),return,end

% get all the display surfaces that match
[topo d.leftDisplaySurfaceNames] = listCaretDir(sprintf('dirname=%s',atlasDir),sprintf('numNodes=%i',d.leftAtlasTopo.numNodes),sprintf('matchStr=LEFT'),'noDisplay=1');
if isempty(d.leftDisplaySurfaceNames),disp(sprintf('(mlrImportCaret) Could not find any display coordinages for LEFT in atlas dir: %s',atlasDir)),return,end
[topo d.rightDisplaySurfaceNames] = listCaretDir(sprintf('dirname=%s',atlasDir),sprintf('numNodes=%i',d.rightAtlasTopo.numNodes),sprintf('matchStr=RIGHT'),'noDisplay=1');
if isempty(d.rightDisplaySurfaceNames),disp(sprintf('(mlrImportCaret) Could not find any display coordinages for RIGHT in atlas dir: %s',atlasDir)),return,end

% get the volume header
[hdr d.volumeFileName] = getNiftiVolumeHeader(surfRelaxDir);

% get how much to shift coords to be in the center
d.coordShift = hdr.dim(2:4)/2;
  
% display the transforms
if verbose
  dispXform('Left hemisphere 711-2B to original xform',d.leftXform);
  dispXform('Right hemisphere 711-2B to original xform',d.rightXform);
end

% Now load left and right coord surfaces
d.leftCoordSurf = loadSurfCaret(fullfile(atlasDir,d.leftAtlasAligned.name),fullfile(atlasDir,d.leftAtlasTopo.name),'xform',d.leftXform,'zeroBased=1','coordShift',d.coordShift);
d.rightCoordSurf = loadSurfCaret(fullfile(atlasDir,d.rightAtlasAligned.name),fullfile(atlasDir,d.rightAtlasTopo.name),'xform',d.rightXform,'zeroBased=1','coordShift',d.coordShift);

% now load display surfaces
d = loadDisplaySurfaces(d,atlasDir);

% compute curvature
d = getCurvature(d,doCurvature,atlasDir);

% check for save directory
if ~isdir('caret'),mkdir('caret');end

% and write out left and right coord surfaces
writeOFF(d.leftCoordSurf,fullfile('caret','leftCoords'));
writeOFF(d.rightCoordSurf,fullfile('caret','rightCoords'));

% write display surfaces
for i = 1:d.displaySurfIndex
  writeOFF(d.displaySurf{i},fullfile('caret',d.displaySurfFileName{i}));
end

% save curvature
if ~isempty(d.leftAtlasCurvature),saveVFF(fullfile('caret','leftAtlasCurvature'),d.leftAtlasCurvature);end
if ~isempty(d.rightAtlasCurvature),saveVFF(fullfile('caret','rightAtlasCurvature'),d.rightAtlasCurvature);end

% and link anatomy file
linkFile(d.volumeFileName,fullfile('caret',getLastDir(d.volumeFileName)));
linkFile(setext(d.volumeFileName,'img'),fullfile('caret',setext(getLastDir(d.volumeFileName),'img')));

%%%%%%%%%%%%%%%%%%%%%%
%%   getCurvature   %%
%%%%%%%%%%%%%%%%%%%%%%
function d = getCurvature(d,doCurvature,atlasDir)

d.leftAtlasCurvature = [];d.rightAtlasCurvature = [];
if doCurvature
  if ~isempty(d.leftFiducialSurfIndex)
    % left curvature filename
    leftCurvatureFile = fullfile(atlasDir,'leftAtlasCurvature.vff');
    % see if we can load it
    if isfile(leftCurvatureFile) && (doCurvature <= 1)
      d.leftAtlasCurvature = loadVFF(leftCurvatureFile);
    else
      % otherwise compute it
      d.leftAtlasCurvature = -calcCurvature(d.displaySurf{d.leftFiducialSurfIndex});
      saveVFF(leftCurvatureFile,d.leftAtlasCurvature);
    end
  end
  if ~isempty(d.rightFiducialSurfIndex)
    % left curvature filename
    rightCurvatureFile = fullfile(atlasDir,'rightAtlasCurvature.vff');
    % see if we can load it
    if isfile(rightCurvatureFile) && (doCurvature <= 1)
      d.rightAtlasCurvature = loadVFF(rightCurvatureFile);
    else
      % otherwise compute it
      d.rightAtlasCurvature = -calcCurvature(d.displaySurf{d.rightFiducialSurfIndex});
      saveVFF(rightCurvatureFile,d.rightAtlasCurvature);
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   loadDisplaySurfaces   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = loadDisplaySurfaces(d,atlasDir)

% names of display surfaces to try to load
surfaceNames = {'fiducial','.inflated','very_Inflated'};
d.displaySurfIndex = 0;d.leftFiducialSurfIndex = [];d.rightFiducialSurfIndex = [];
for i = 1:length(surfaceNames)
  % try to load left
  surf = myLoadSurfCaret(d,atlasDir,'left',surfaceNames{i});
  % save it if we got it.
  if ~isempty(surf)
    d.displaySurfIndex = d.displaySurfIndex+1;
    d.displaySurf{d.displaySurfIndex} = surf;
    % create name for file
    d.displaySurfFileName{d.displaySurfIndex} = fixBadChars(surfaceNames{i},{{'.',''},{'_',''}});
    d.displaySurfFileName{d.displaySurfIndex} = sprintf('leftAtlas%s',capitalize(d.displaySurfFileName{d.displaySurfIndex}));
    % if this is the first, fiducial one, then remember it
    if i == 1,d.leftFiducialSurfIndex = d.displaySurfIndex;end
  end
  % try to load right
  surf = myLoadSurfCaret(d,atlasDir,'right',surfaceNames{i});
  % save it if we got it.
  if ~isempty(surf)
    d.displaySurfIndex = d.displaySurfIndex+1;
    d.displaySurf{d.displaySurfIndex} = surf;
    % create name for file
    d.displaySurfFileName{d.displaySurfIndex} = fixBadChars(surfaceNames{i},{{'.',''},{'_',''}});
    d.displaySurfFileName{d.displaySurfIndex} = sprintf('rightAtlas%s',capitalize(d.displaySurfFileName{d.displaySurfIndex}));
    % if this is the first, fiducial one, then remember it
    if i == 1,d.rightFiducialSurfIndex = d.displaySurfIndex;end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%   myLoadSurfCaret   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function surf = myLoadSurfCaret(d,atlasDir,side,type)

surf = [];
% get the surface name
surfaceName = dirfind(d.(sprintf('%sDisplaySurfaceNames',side)),type);
% check if we have one
if isempty(surfaceName)
  disp(sprintf('(mlrImportCaret:myLoadSurfCaret) Could not find %s %s atlas coord file in %s',side,type,atlasDir));
  return
end
% load it 
surf = loadSurfCaret(fullfile(atlasDir,surfaceName),fullfile(atlasDir,d.(sprintf('%sAtlasTopo',side)).name));
%%%%%%%%%%%%%%%%%%%%%%
%%   getAtlasTopo   %%
%%%%%%%%%%%%%%%%%%%%%%
function [leftAtlasTopo rightAtlasTopo] = getAtlasTopo(atlasDir)

leftAtlasTopo = listCaretDir(sprintf('dirname=%s',atlasDir),'matchStr=LEFT','noDisplay=1');
if isempty(leftAtlasTopo)
  disp(sprintf('(mlrImportCaret) Could not find any atlas topo files that match LEFT in %s',atlasDir));
  return
end
rightAtlasTopo = listCaretDir(sprintf('dirname=%s',atlasDir),'matchStr=RIGHT','noDisplay=1');
if isempty(rightAtlasTopo)
  disp(sprintf('(mlrImportCaret) Could not find any atlas topo files that match RIGHT in %s',atlasDir));
  return
end

% ask user to select topo files
paramsInfo{1} = {'leftTopo',dir2cell(leftAtlasTopo),'type=popupmenu','Choose the topo file for the left hemisphere'};
paramsInfo{2} = {'rightTopo',dir2cell(rightAtlasTopo),'type=popupmenu','Choose the topo file for the left hemisphere'};
params = mrParamsDialog(paramsInfo,'Select atlas topo files');
if isempty(params),return,end

% get the topo file info
for i = 1:length(leftAtlasTopo)
  if isequal(leftAtlasTopo(i).name,params.leftTopo)
    leftAtlasTopo = leftAtlasTopo(i);
    break;
  end
end
% get the topo file info
for i = 1:length(rightAtlasTopo)
  if isequal(rightAtlasTopo(i).name,params.rightTopo)
    rightAtlasTopo = rightAtlasTopo(i);
    break;
  end
end

%%%%%%%%%%%%%%%%%%%
%    dispXform    %
%%%%%%%%%%%%%%%%%%%
function dispXform(name,xform)

disp(sprintf('(importCaret:dispXform) %s',name));
for i = 1:size(xform,1)
  disp(sprintf('%s',mynum2str(xform(i,:),'sigfigs=4')));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    computeTransformBetweenSurfaces    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xform = computeTransformBetweenSurfaces(fromSurfaceFilename,toSurfaceFilename)

xform = [];
% load the surfaces
fromSurface = openCaretFile(fromSurfaceFilename);
if isempty(fromSurface),return,end
toSurface = openCaretFile(toSurfaceFilename);
if isempty(toSurface),return,end

% check matching nodes
if fromSurface.num_nodes ~= toSurface.num_nodes
  disp(sprintf('(importCaret:computeTransformBetweenSurfaces) Nodes must match between 711-2B surface (%s:%i) and original surface (%s:%i)',fromSurfaceFilename,fromSurface.num_nodes,toSurfaceFilename,toSurface.num_nodes));
  return
end

% now compute the best matching xform between these two surfaces
toVertices = toSurface.data';toVertices(4,:) = 1;
fromVertices = fromSurface.data';fromVertices(4,:) = 1;
xform = toVertices*pinv(fromVertices);
xform(4,:) = [0 0 0 1];

%%%%%%%%%%%%%%%%%%%%%%%
%    getCoordFiles    %
%%%%%%%%%%%%%%%%%%%%%%%
function [filenames stemName] = getCoordFiles(caretFileDir,matchStr)

filenames = [];
[topoFile coordFile] = listCaretDir(sprintf('dirname=%s',caretFileDir),sprintf('matchStr=%s.coord',matchStr),'noDisplay=1');
if isempty(coordFile)
  disp(sprintf('(mlrImportCaret) Could not find %s.coord file',matchStr));
  return
end
if length(coordFile) ~= 1
  disp(sprintf('(mlrImportCaret) Found %i files that match %s.coord file, but should have found only 1',matchStr));
  return
end  

% now go through and construct rest of files and make sure they are there
matchLoc = findstr(coordFile.name,matchStr);
stemName = sprintf('%s%s',coordFile.name(1:matchLoc-1),matchStr);

% construct the name of the surface which has been transformed into
% the Caret 711-2B space.
caretName = sprintf('%s_711-2B.coord',stemName);
caretName = fullfile(caretFileDir,caretName);
if ~isfile(caretName)
  disp(sprintf('(mlrImportCaret) Could not find file %s',caretName));
  return
end

filenames{1} = fullfile(caretFileDir,coordFile.name);
filenames{2} = caretName;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    checkAllCaretFiles    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkAllCaretFiles(surfRelaxDir,caretFileDir)

% check directory
disppercent(-inf,sprintf('(mlrImportCaret) Checking %s directory for topo and coord files',caretFileDir));
[topoFiles coordFiles] = listCaretDir(sprintf('dirname=%s',caretFileDir),'noDisplay=1');
disppercent(inf);

% return if files not found
if isempty(topoFiles)
  disp(sprintf('(mlrImportCaret) Could not find any topo files in %s. Have you run mypreborder/mypostborder?',caretFileDir));
  return
end
% return if files not found
if isempty(coordFiles)
  disp(sprintf('(mlrImportCaret) Could not find any coord files in %s. Have you run mypreborder/mypostborder?',caretFileDir));
  return
end

% get topoFileInfo
for i = 1:length(topoFiles)
  topoFileNodes{i} = topoFiles(i).numNodes;
  topoFileComment{i} = topoFiles(i).comment;
end

% get coordFileInfo
for i = 1:length(coordFiles)
  coordFileNodes{i} = coordFiles(i).numNodes;
  coordFileComment{i} = coordFiles(i).comment;
end

% get nifti volume files
niftiVolumes = dir(fullfile(surfRelaxDir,'*.hdr'));
niftiVolumes = cat(1,niftiVolumes,dir(fullfile(surfRelaxDir,'*.nii')));

% check that we got something
if isempty(niftiVolumes)
  disp(sprintf('(mlrImportCaret) Could not find any nifti volume files in %s',surfRelaxDir));
  return
end

paramsInfo{1} = {'topoFileNum',1,'incdec=[-1 1]',sprintf('minmax=[1 %i]',length(topoFiles)),'Topo files contain the topology, i.e. the list of triangles that tile the surface'};
paramsInfo{end+1} = {'topoFileName',dir2cell(topoFiles),'type=String','group=topoFileNum','editable=0','Topo filename'};
paramsInfo{end+1} = {'topoFileNodes',topoFileNodes,'group=topoFileNum','type=numeric','editable=0','Number of nodes in topo file. Must match coord file.'};
paramsInfo{end+1} = {'topoFileComment',topoFileComment,'group=topoFileNum','type=string','editable=0','Comment for the topo files'};
paramsInfo{end+1} = {'coordFileNum',1,'incdec=[-1 1]',sprintf('minmax=[1 %i]',length(coordFiles)),'Coord files contain the coordinates of each vertex in the surface'};
paramsInfo{end+1} = {'coordFileName',dir2cell(coordFiles),'type=String','group=coordFileNum','editable=0','Name of coord file'};
paramsInfo{end+1} = {'coordFileNodes',coordFileNodes,'type=numeric','group=coordFileNum','editable=0','Number of nodes in coord file. Must match topo file.'};
paramsInfo{end+1} = {'coordFileComment',coordFileComment,'type=String','group=coordFileNum','editable=0','Comment for the coord file'};
paramsInfo{end+1} = {'niftiVolume',dir2cell(niftiVolumes),'The nifti volume should be the canonical volume from which this surface was made'};

% bring up dialog box
params = mrParamsDialog(paramsInfo,'Select topo and coord file to display');
if isempty(params),return,end

% check matching nodes
if topoFiles(params.topoFileNum).numNodes ~= coordFiles(params.coordFileNum).numNodes
  disp(sprintf('(mlrImportCaret) Num nodes between coord (%i) and topo (%i) do not match',coordFiles(params.coordFileNum).numNodes,topoFiles(params.topoFileNum).numNodes));
  return
end

% get the files
topoFileName = fullfile(caretFileDir,params.topoFileName{params.topoFileNum});
coordFileName = fullfile(surfRelaxDir,params.coordFileName{params.coordFileNum});
niftiFileName = fullfile(surfRelaxDir,params.niftiVolume);

% load the surface
surf = loadSurfCaret(coordFileName,topoFileName);
keyboard

%%%%%%%%%%%%%%%%%
%    dirfind   %%
%%%%%%%%%%%%%%%%%
function filename = dirfind(d,matchstr)

filename = '';
for i = 1:length(d)
  if ~isempty(strfind(lower(d(i).name),lower(matchstr)))
    filename = d(i).name;
    return
  end
end

    
%%%%%%%%%%%%%%%%%%
%    dir2cell    %
%%%%%%%%%%%%%%%%%%
function c = dir2cell(d)
		    
for i = 1:length(d)
  c{i} = d(i).name;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getNiftiVolumeHeader    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hdr volumeFileName] = getNiftiVolumeHeader(surfRelaxDir)

hdr = [];
volumeFileName = [];
% get nifti volume files
niftiVolumes = dir(fullfile(surfRelaxDir,'*.hdr'));
niftiVolumes = cat(1,niftiVolumes,dir(fullfile(surfRelaxDir,'*.nii')));
% check that we got something
if isempty(niftiVolumes)
  disp(sprintf('(mlrImportCaret:getNiftiVolumeHeader) Could not find any nifti volume files in %s',surfRelaxDir));
  return
end
if length(niftiVolumes) > 1
  paramsInfo{1} = {'niftiVolume',dir2cell(niftiVolumes),'The nifti volume should be the canonical volume from which this surface was made'};

  % bring up dialog box
  params = mrParamsDialog(paramsInfo,'Select topo and coord file to display');
  if isempty(params),return,end
  volumeFileName = fullfile(surfRelaxDir,params.niftiVolume);
  hdr = mlrImageReadNiftiHeader(volumeFileName);
else
  volumeFileName = fullfile(surfRelaxDir,niftiVolumes.name);
  hdr = mlrImageReadNiftiHeader(volumeFileName);
end


%%%%%%%%%%%%%%%%%%%
%%   checkDirs   %%
%%%%%%%%%%%%%%%%%%%
function isok = checkDirs(varargin)

isok = 0;
for i = 1:nargin
  if ~isdir(varargin{i})
    disp(sprintf('(mlrImportCaret) Could not find %s directory',varargin{i}));
    return
  end
end
isok = 1;

%%%%%%%%%%%%%%%%%%%%
%%   capitalize   %%
%%%%%%%%%%%%%%%%%%%%
function name = capitalize(name)

name(1) = upper(name(1));
