% importFlatOFF.m
%
%        $Id$	
%      usage: importFlatOFF(v,'flatPatch.off')
%         by: eli merriam
%       date: 09/11/07
%    purpose: load a surfRelax flatpatch, and everything that goes with it
%  
%    Three ways to use this code:
function base = importFlatOFF(flatFile)

% check arguments
if ~any(nargin == [0 1 2])
  help importFlatOFF
  return;
end
base = [];
if ~exist('mrParamsDialog');
  disp(sprintf('(importFlatOFF) You must have mrTools in your path to run this'));
  return;
end

% init arguments
if nargin == 0
  % Open dialog box to have user choose the file
  startPathStr = mrGetPref('volumeDirectory');
  filterspec = {'*.off','SurfRelax off file';'*lat.off','SurfRelax off flat file';'*.*','All files'};
  title = 'Choose flat OFF file';
  flatFile = mlrGetPathStrDialog(startPathStr,title,filterspec,'off');
  flatFile = cellArray(flatFile);
  flatFile = flatFile{1};
end

if isstr(flatFile);
  % check to see if we are passed in a file name
  if ~isfile(flatFile);
    disp(sprintf('(importFlatOFF) %s does not exist', flatFile));
    return;
  end
  % load with mrFlatViewer
  params = mrFlatViewer(flatFile);
  if isempty(params),return,end

  % mrFlatViewer  will create the gFlatViewer global structure
  % we can then use that for a lot of the further processing
  global gFlatViewer
end


% could be passed in a params structure, rather than a flat file
% this is the case when importFlatOFF is called from makeFlat
if isstruct(flatFile);
  params = flatFile;
  if isfile(fullfile(params.path, params.flatFileName));
    flatFile = loadSurfOFF(fullfile(params.path, params.flatFileName));
    surf.inner = loadSurfOFF(fullfile(params.path, params.innerCoordsFileName));
    surf.outer = loadSurfOFF(fullfile(params.path, params.outerCoordsFileName));
    surf.curv = loadVFF(fullfile(params.path, params.curvFileName));
    anat.hdr = mlrImageReadNiftiHeader(fullfile(params.path, params.anatFileName));
    % now do necessary world2array xformation here
    surf.outer = xformSurfaceWorld2Array(surf.outer,anat.hdr);
    surf.inner = xformSurfaceWorld2Array(surf.inner,anat.hdr);
  end
end

% just crap out if the params still doesn't exist
if ieNotDefined('params');
  return;
end


% get two additional parameters if they were not passed in before
if ~any(isfield(params, {'threshold', 'flatRes'}));
  paramsInfo = {};
  paramsInfo{end+1} = {'threshold', 1, 'type=checkbox', 'Whether or not to threshold the flat patch'};
  paramsInfo{end+1} = {'flipFlag', 0, 'type=checkbox', 'Some patches come out flipped.  Use this option to correct this problem'};
  paramsInfo{end+1} = {'flatRes', 2, 'incdec=[-1 1]', 'Factore by which the resolution of the flat patch is increased'};
  paramsInfo{end+1} = {'flatBaseName', stripext(getLastDir(params.flatFileName)), 'Name of the flat base anatomy'};
  flatParams = mrParamsDialog(paramsInfo,'Flat patch parameters');
  % check for cancel
  if isempty(flatParams);
    return;
  end
else
  flatParams.threshold = params.threshold;
  flatParams.flatRes = params.flatRes;
  flatParams.flipFlag = 0;
  flatParams.flatBaseName = stripext(getLastDir(params.flatFileName));
end


% check to see if we got here from the flatViewer
% we need to translate a few variable names if we did
if ~ieNotDefined('gFlatViewer');
  flat.whichInx  = gFlatViewer.flat.patch2parent(:,2); 
  flat.locsFlat  = gFlatViewer.flat.vtcs;
  flat.curvature = gFlatViewer.curv(flat.whichInx);
  flat.locsInner = gFlatViewer.surfaces.inner.vtcs(flat.whichInx,:);
  flat.locsOuter = gFlatViewer.surfaces.outer.vtcs(flat.whichInx,:);
  flat.hdr       = gFlatViewer.anat.hdr;
elseif ~isempty(flatFile);
  % or maybe we got here from loading a flatFile
  flat.whichInx  = flatFile.patch2parent(:,2); 
  flat.locsFlat  = flatFile.vtcs;
  flat.locsInner = surf.inner.vtcs(flat.whichInx,:);
  flat.locsOuter = surf.outer.vtcs(flat.whichInx,:);
  flat.curvature = surf.curv(flat.whichInx);
  flat.hdr       = anat.hdr;
else
  disp(sprintf('(importFlatOFF) cannot find parameters needed to make base anatomy'));
end  


% this X-Y swaping only changes the orientation of the image
% isn't a crucial step
flat.locsFlat = [flat.locsFlat(:,2) flat.locsFlat(:,1) flat.locsFlat(:,3)];

flat.minLocsFlat = min(flat.locsFlat);
flat.locsFlat(:,1) = flat.locsFlat(:,1) - flat.minLocsFlat(1) + 1;
flat.locsFlat(:,2) = flat.locsFlat(:,2) - flat.minLocsFlat(2) + 1;

imSize = round(max(flat.locsFlat));

if flatParams.flipFlag == 1
  disp(sprintf('(importFlatOFF) flipFlag was set to 1, so X-Y flipping flat patch'))
  x = flat.locsFlat(:,2);
  y = flat.locsFlat(:,1);
else
  x = flat.locsFlat(:,1);
  y = flat.locsFlat(:,2);
end

xi = [1:(1/flatParams.flatRes):imSize(1)];
yi = [1:(1/flatParams.flatRes):imSize(2)]';

% why is this flip required???  -epm
yi = flipud(yi);

warning off
flat.map = griddata(x,y,flat.curvature,xi,yi,'linear');

% grid the 3d coords
for i=1:3
  flat.baseCoordsInner(:,:,i) =  griddata(x,y, flat.locsInner(:,i), xi, yi,'linear');
  flat.baseCoordsOuter(:,:,i) =  griddata(x,y, flat.locsOuter(:,i), xi, yi,'linear');
end

warning on
% mask out out-of-brain coords
flat.baseCoordsInner(isnan(flat.map)) = 0;
flat.baseCoordsOuter(isnan(flat.map)) = 0;

% get rid of any NaN's
flat.baseCoordsInner(~isfinite(flat.baseCoordsInner)) = 0;
flat.baseCoordsOuter(~isfinite(flat.baseCoordsOuter)) = 0;

% base.map = mrUpSample(base.map);

% now blur image
flat.blurMap(:,:) = blur(flat.map(:,:));
% threshold image
flat.thresholdMap(:,:) = (flat.blurMap(:,:)>nanmedian(flat.blurMap(:)))*0.5+0.25;
% flat.medianVal = nanmedian(flat.blurMap(:));
% flat.blurMap(flat.blurMap<flat.medianVal) = -1;
% flat.blurMap(flat.blurMap>flat.medianVal) = 1;
flat.thresholdMap = blur(flat.thresholdMap);
flat.thresholdMap(isnan(flat.map)) = 0;

% now load the anatomy file, so that we can get the vol2mag/vol2tal fields
v = newView;
%if extension is hdr, change to img
[~,~,anatFileExtension] = fileparts(params.anatFileName);
if strcmp(anatFileExtension,'hdr')
  anatfileName = setext(params.anatFileName,'img');
else
  anatfileName = params.anatFileName;
end
v = loadAnat(v,anatfileName,params.path);
b = viewGet(v,'base');
deleteView(v);

% now generate a base structure
clear base;
base.hdr = flat.hdr;
base.name = flatParams.flatBaseName;
base.vol2mag = b.vol2mag;
base.vol2tal = b.vol2tal;

% Extract permutation matrix to keep track of slice orientation.
base.permutationMatrix = getPermutationMatrix(base.hdr);

% set base parameters
base.coordMap.path = params.path;
base.coordMap.flatFileName = params.flatFileName;
base.coordMap.innerCoordsFileName = params.innerCoordsFileName;
base.coordMap.outerCoordsFileName = params.outerCoordsFileName;
base.coordMap.curvFileName = params.curvFileName;
base.coordMap.anatFileName = params.anatFileName;

% load all the flat maps into the base. We
% need to make all the flat images have
% the same width and height.
if flatParams.threshold
  base.data(:,:,1) = flat.thresholdMap;
  base.range = [0 1];
  base.clip = [0 1];
else
  flat.map(flat.map>1) = 1;
  flat.map(flat.map<-1) = -1;
  flat.map = 32*(flat.map+1);
  base.data(:,:,1) = flat.map;
end    

% start of coords as inner coords 
base.coordMap.coords(:,:,1,1) = flat.baseCoordsInner(:,:,1);
base.coordMap.coords(:,:,1,2) = flat.baseCoordsInner(:,:,2);
base.coordMap.coords(:,:,1,3) = flat.baseCoordsInner(:,:,3);
% setup inner coords
base.coordMap.innerCoords(:,:,1,1) = flat.baseCoordsInner(:,:,1);
base.coordMap.innerCoords(:,:,1,2) = flat.baseCoordsInner(:,:,2);
base.coordMap.innerCoords(:,:,1,3) = flat.baseCoordsInner(:,:,3);
% setup outer coords
base.coordMap.outerCoords(:,:,1,1) = flat.baseCoordsOuter(:,:,1);
base.coordMap.outerCoords(:,:,1,2) = flat.baseCoordsOuter(:,:,2);
base.coordMap.outerCoords(:,:,1,3) = flat.baseCoordsOuter(:,:,3);

% set dimensions
base.coordMap.dims = flat.hdr.dim([2 3 4])';

% set base
base.range = [min(min(base.data)) max(max(base.data))];
base.clip = [0 1.5];
base.type = 1;

clear global gFlatViewer;

