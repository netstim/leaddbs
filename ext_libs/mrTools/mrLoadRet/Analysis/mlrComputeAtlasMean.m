% mlrComputeAtlasMean.m
%
%        $Id$ 
%      usage: mlrComputeAtlasMean(sessionPath,groupNum,scanNum,analysisName,overlayName,baseAnatomy)
%         by: justin gardner
%       date: 12/11/09
%    purpose: To compute a mean overlay across subjects. The mean overlay is taken for a particular
%             baseAnatomy. So, if you want to take a mean across subjects, the baseAnatomy should
%             be an atlas brain like one imported from Caret (see mlrImportCaret). A new overlay
%             will be saved in which one of the sessions as mrDispOverlayAnal/meanAnalysis with
%             an overlay called meanOveraly. Note that the analysis and baseAnatomy that you are
%             using has to be saved in each one of the sessions (i.e. the base anatomy has to be in
%             the Anatomy directory -- use File/Base Anatomy/Save) and the analysis has to be
%             in an appropriate directory under the Group (use File/Analysis/Save if it is not).
%
%             e.g. 
%             mlrComputeAtlasMean('leftAtlasVeryInflated',{'S00320090717','S00920090717'},'groupNum=Concatenation','scanNum=1','analysisName=erAnal','overlayName=r2');
%
%             Note that you can specify particular parameters specifically for each session, or
%             just have all of the sessions share the same parameters, as above. For example,
%             if you wanted to use the second scan instead of the first scan for sesson 2, above:
%
%             mlrComputeAtlasMean('leftAtlasVeryInflated',{'S00320090717','S00920090717'},'groupNum=Concatenation','scanNum=[1 2]','analysisName=erAnal','overlayName=r2');
%
%             You can also do more than one base anatomy at a time:
% 
%             mlrComputeAtlasMean({'leftAtlasVeryInflated','rightAtlasVeryInflated'},{'S00320090717','S00920090717'},'groupNum=Concatenation','scanNum=1','analysisName=erAnal','overlayName=r2');
%
%             You can also compute a "mean roi"
%
%             mlrComputeAtlasMean('leftAtlasVeryInflated',{'S00320090717','S00920090717'},'roiName=fef');
%             You have the option of computing the intersect (default) of ROIs or the union:
%             mlrComputeAtlasMean('leftAtlasVeryInflated',{'S00320090717','S00920090717'},'roiName=fef','roiMergeFun=union');
function retval = mlrComputeAtlasMean(baseAnatomies, sessionPath, varargin)

% check arguments
if nargin < 3
  help mlrComputeAtlasMean
  return
end

groupNum = [];
scanNum = [];
analysisName = [];
overlayName = [];
roiName = [];
roiMergeFun = [];
getArgs(varargin,{'groupNum=[]','scanNum=[]','analysisName=[]','overlayName=[]','roiName=[]','roiMergeFun=intersect'});

% parse input arguments
sessions = parseArguments(baseAnatomies,sessionPath,groupNum,scanNum,analysisName,overlayName,roiName);
if isempty(sessions),return,end

% do the overlays
if sessions.doOverlay
  % go through each session and load the overlays
  overlays = loadOverlays(sessions);
  if isempty(overlays),return,end

  % compute the average overlays
  meanOverlays = computeMeanOverlays(sessions,overlays);

  % save the overlays back
  saveMeanOverlays(sessions,meanOverlays);
end

% do the rois
if sessions.doROI
  % load the ROIs
  ROIs = loadROIs(sessions);
  if isempty(ROIs),return,end

  % compute the average roi
  meanROIs = computeMeanROIs(sessions,ROIs,roiMergeFun);
  if isempty(meanROIs),return,end
  
  % save the overlays back
  saveMeanROIs(sessions,meanROIs,roiMergeFun);
end


%%%%%%%%%%%%%%%%%%%%%%
%%   saveMeanROIs   %%
%%%%%%%%%%%%%%%%%%%%%%
function saveMeanROIs(sessions,meanROIs,roiMergeFun)

for iSession = 1:sessions.nSessions
  % open a view on to the session
  thisPwd = pwd;
  cd(sessions.sessionPath{iSession});
  v = newView;
  cd(thisPwd);
  % now project the overlay back into the scan coordinates
  % load the anatomies 
  for iBase = 1:sessions.nBases
    % load the anatomy
    v = loadAnat(v,sessions.baseAnatomies{iBase});
    % get the base
    b = viewGet(v,'baseCoordMap');
    %we'll take the coordinates of the middle of whatever range of cortical depth is currenlty selected
    corticalSlice = ceil(mean(viewGet(v,'corticalDepth'))*size(baseCoordMap.coords,5));
    coords = squeeze(b.coords(:,:,:,:,corticalSlice));
    
    % and put the coordinates of those vertices into our search ROI
    thisROI.coords = coords(meanROIs{iBase}.vertices,:)';
    thisROI.coords(4,:) = 1;

    % set other parts of the search ROI
    thisROI.name = sprintf('%sAtlas%s%s',sessions.roiName{iSession},upper(roiMergeFun(1)),lower(roiMergeFun(2:end)));
    thisROI.voxelSize=viewGet(v,'baseVoxelSize');
    thisROI.xform=viewGet(v,'basesform');
    [tf thisROI] = isroi(thisROI);
    
    % and install the roi
    v = viewSet(v,'newROI',thisROI,1);
    saveROI(v,thisROI.name);
  end
  % close and delete
  deleteView(v);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%   computeMeanROIs   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function meanROI = computeMeanROIs(sessions,ROIs,roiMergeFun)

meanROI = {};
for iBase = 1:sessions.nBases
  for iSession = 1:sessions.nSessions
    % sum up overlays across all sessions
    if iSession == 1
      meanROI{iBase}.vertices = ROIs.r{iSession,iBase}.vertices;
    else
      if strcmp(roiMergeFun,'intersect')
	meanROI{iBase}.vertices = intersect(meanROI{iBase}.vertices,ROIs.r{iSession,iBase}.vertices);
      elseif strcmp(roiMergeFun,'union')
	meanROI{iBase}.vertices = union(meanROI{iBase}.vertices,ROIs.r{iSession,iBase}.vertices);
      else
	disp(sprintf('(mlrComputeAtlasMean:computeMeanROI) Unknown roiMergeFun %s',roiMergeFun));
	meanROI = {};
	return
      end
    end
  end
  % and display
  dispROI(ROIs.b{iSession,iBase},cellcat(ROIs.r,meanROI{iBase}),sprintf('Mean on %s',sessions.baseAnatomies{iBase}));
end


%%%%%%%%%%%%%%%%%%
%    loadROIs    %
%%%%%%%%%%%%%%%%%%
function ROIs = loadROIs(sessions)

ROIs = [];

for iSession = 1:sessions.nSessions
  mrQuit([]);
  % first make sure the directory exists
  if ~myisdir(sessions.sessionPath{iSession})
    disp(sprintf('(mlrComputeAtlasMean:loadROIs) Could not find directory %s',sessions.sessionPath{iSession}));
    return
  end
  % next make sure it has a session
  if ~isfile(fullfile(sessions.sessionPath{iSession},'mrSession.mat'))
    disp(sprintf('(mlrComputeAtlasMean:loadROIs) Cound not find mrSession.mat in %s',sessions.sessionPath{iSession}));
    return
  end
  % create a view for each session
  thisPwd = pwd;
  cd(sessions.sessionPath{iSession});
  v = newView;
  cd(thisPwd);
  % load the roi
  v = loadROI(v,sessions.roiName{iSession});
  if viewGet(v,'numROIs') < 1,return,end
  % load the anatomies 
  for iBase = 1:sessions.nBases
    % load the anatomy
    v = loadAnat(v,sessions.baseAnatomies{iBase});
    % get the base
    b{iSession,iBase} = viewGet(v,'baseSurface');
    if isempty(b(iSession,iBase)), disp(sprintf('(mlrComputeAtlasMean:loadROIs) Could not load base'));return,end
    % get the default overlay
    [b{iSession,iBase}.overlayImage base thisROI] = refreshMLRDisplay(viewGet(v,'viewNum'));
    % get the display vertices
    roi{iSession,iBase}.vertices = thisROI{1}.vertices;
    roi{iSession,iBase}.overlayImage = thisROI{1}.overlayImage;
    % display in a figure
    dispROI(b{iSession,iBase},roi{iSession,iBase},makeSessionName(v,sessions,iSession));
  end
  % delete the view
  deleteView(v);
end

% check that all the bases match
for iBase = 1:sessions.nBases
  for iSession = 2:sessions.nSessions
    if size(b{1,iBase}.vtcs,1) ~= size(b{iSession,iBase}.vtcs,1) 
      disp(sprintf('(mlrComputeAtlasMean:loadOverlays) Base %s does not match number of vertices for %s',sessions.baseAnatomies{iBase},makeSessionName([],sessions,iSession)));
      return
    end
  end
end
% set return fields
ROIs.b = b;
ROIs.r = roi;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   saveMeanOverlays   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveMeanOverlays(sessions,meanOverlays)

for iSession = 1:sessions.nSessions
  % open a view on to the session
  thisPwd = pwd;
  cd(sessions.sessionPath{iSession});
  v = newView;
  cd(thisPwd);
  % set the group
  v = viewSet(v,'curGroup',sessions.groupNum{iSession});
  % set the scan number
  v = viewSet(v,'curScan',sessions.scanNum{iSession});
  % get scan dims
  scanDims = viewGet(v,'scanDims');
  % create an empty overlay
  overlay = nan(scanDims);
  % now project the overlay back into the scan coordinates
  % load the anatomies 
  for iBase = 1:sessions.nBases
    % load the anatomy
    v = loadAnat(v,sessions.baseAnatomies{iBase});
    % get the base
    b = viewGet(v,'baseCoordMap');
    % get the transformation from base2scan coordinates
    base2scan = viewGet(v,'base2scan');
    % get the coordinates of each vertex
    %we'll take the coordinates of the middle of whatever range of cortical depth is currenlty selected
    corticalSlice = ceil(mean(viewGet(v,'corticalDepth'))*size(baseCoordMap.coords,5));
    baseCoords = squeeze(b.coords(:,:,:,:,corticalSlice));
    baseCoords(:,4) = 1;
    % convert into scan coords
    scanCoords = round(base2scan*baseCoords');
    scanLinearCoords = mrSub2ind(scanDims,scanCoords(1,:),scanCoords(2,:),scanCoords(3,:));
    % get the overlay points that lay within the image
    overlayIm = meanOverlays{iBase}.overlayIm;
    % get points outside of scan and remove
    goodPoints = find(~isnan(scanLinearCoords));
    scanLinearCoords = scanLinearCoords(goodPoints);
    overlayIm = overlayIm(goodPoints);
    % and put the overlay into the scan coordinate iamge
    overlay(scanLinearCoords) = overlayIm;
  end
  % and install the overlay into a custom analysis
  mrDispOverlay(overlay,sessions.scanNum{iSession},viewGet(v,'groupNum',sessions.groupNum{iSession}),[],'overlayName=meanOverlay','saveName=meanAnalysis');
  % close and delete
  deleteView(v);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    computeMeanOverlay    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function meanOverlays = computeMeanOverlays(sessions,overlays)

% NOTE: This averaging won't work for amp/phase maps like the corAnal which have
% to be converted back to a complex number and averaged.

for iBase = 1:sessions.nBases
  for iSession = 1:sessions.nSessions
    % sum up overlays across all sessions
    if iSession == 1
      meanOverlays{iBase}.overlayIm = overlays.o(iSession,iBase).overlayIm;
    else
      meanOverlays{iBase}.overlayIm = meanOverlays{iBase}.overlayIm + overlays.o(iSession,iBase).overlayIm;
    end
  end
  % divide by n
  meanOverlays{iBase}.overlayIm = meanOverlays{iBase}.overlayIm/sessions.nSessions;
  % convert to RGB
  %clim = [min(meanOverlays{iBase}.overlayIm) max(meanOverlays{iBase}.overlayIm)];
  clim = overlays.o(iSession,iBase).range;
  meanOverlays{iBase}.overlayRGB = rescale2rgb(meanOverlays{iBase}.overlayIm,overlays.o(iSession,iBase).cmap,clim);
  % and display
  dispOverlay(overlays.b(iSession,iBase),meanOverlays{iBase}.overlayRGB,sprintf('Mean on %s',sessions.baseAnatomies{iBase}));
end


%%%%%%%%%%%%%%%%%%%%%%
%    loadOverlays    %
%%%%%%%%%%%%%%%%%%%%%%
function overlays = loadOverlays(sessions)

overlays = [];

for iSession = 1:sessions.nSessions
  mrQuit([]);
  % first make sure the directory exists
  if ~myisdir(sessions.sessionPath{iSession})
    disp(sprintf('(mlrComputeAtlasMean:loadOverlays) Could not find directory %s',sessions.sessionPath{iSession}));
    return
  end
  % next make sure it has a session
  if ~isfile(fullfile(sessions.sessionPath{iSession},'mrSession.mat'))
    disp(sprintf('(mlrComputeAtlasMean:loadOverlays) Cound not find mrSession.mat in %s',sessions.sessionPath{iSession}));
    return
  end
  % create a view for each session
  thisPwd = pwd;
  cd(sessions.sessionPath{iSession});
  v = newView;
  cd(thisPwd);
  % set the group
  v = viewSet(v,'curGroup',sessions.groupNum{iSession});
  % set the scan number
  v = viewSet(v,'curScan',sessions.scanNum{iSession});
  % load the analysis
  v = loadAnalysis(v,sessions.analysisName{iSession});
  if viewGet(v,'numAnalyses') < 1,return,end
  % load the overlay, just to make sure it is there
  v = viewSet(v,'curOverlay',sessions.overlayName{iSession});
  o{iSession} = viewGet(v,'overlay');
  if isempty(o{iSession}),return,end
  % load the anatomies 
  for iBase = 1:sessions.nBases
    % load the anatomy
    v = loadAnat(v,sessions.baseAnatomies{iBase});
    % get the base
    b(iSession,iBase) = viewGet(v,'baseSurface');
    if isempty(b(iSession,iBase)), disp(sprintf('(mlrComputeAtlasMean:loadOverlays) Could not load base'));return,end
    % compute the overlay
    [overlayImage base roi overlay(iSession,iBase)] = refreshMLRDisplay(viewGet(v,'viewNum'));
    if isempty(overlay(iSession,iBase)),return,end
    % display in a figure
    dispOverlay(b(iSession,iBase),overlay(iSession,iBase).RGB,makeSessionName(v,sessions,iSession));
  end
  % delete the view
  deleteView(v);
end

% check that all the bases match
for iBase = 1:sessions.nBases
  for iSession = 2:sessions.nSessions
    if size(b(1,iBase).vtcs,1) ~= size(b(iSession,iBase).vtcs,1) 
      disp(sprintf('(mlrComputeAtlasMean:loadOverlays) Base %s does not match number of vertices for %s',sessions.baseAnatomies{iBase},makeSessionName([],sessions,iSession)));
      return
    end
  end
end
% set return fields
overlays.b = b;
overlays.o = overlay;

%%%%%%%%%%%%%%%%%%%%%%%%%
%    makeSessionName    %
%%%%%%%%%%%%%%%%%%%%%%%%%
function s = makeSessionName(v,sessions,i)

if sessions.doOverlay
  if isempty(v)
    s = sprintf('%s: %s:%i %s:%s',sessions.sessionPath{i},sessions.groupNum{i},sessions.scanNum{i},sessions.analysisName{i},sessions.overlayName{i});
  else
    s = sprintf('%s: %s:%i %s:%s',sessions.sessionPath{i},viewGet(v,'groupName',sessions.groupNum{i}),sessions.scanNum{i},sessions.analysisName{i},viewGet(v,'overlayName',sessions.overlayName{i}));
  end
else
  s = sprintf('%s',sessions.sessionPath{i});
end
%%%%%%%%%%%%%%%%%%%%%
%    dispOverlay    %
%%%%%%%%%%%%%%%%%%%%%
function dispOverlay(b,overlay,titleStr)

f = mlrSmartfig(fixBadChars(sprintf('mlrComputeAtlasMean:%s',titleStr)),'reuse');clf
patch('vertices', b.vtcs, 'faces', b.tris,'FaceVertexCData', squeeze(overlay),'facecolor','interp','edgecolor','none');
axis equal;
axis off;
camva(4);
rotate3d on;
title(titleStr);
set(gca(f),'XDir','reverse');
set(gca(f),'YDir','normal');
set(gca(f),'ZDir','normal');
drawnow

%%%%%%%%%%%%%%%%%
%    dispROI    %
%%%%%%%%%%%%%%%%%
function dispROI(b,r,titleStr)

f = mlrSmartfig(fixBadChars(sprintf('mlrComputeAtlasMean:%s',titleStr)),'reuse');clf

% draw surface
patch('vertices', b.vtcs, 'faces', b.tris,'FaceVertexCData', squeeze(b.overlayImage),'facecolor','interp','edgecolor','none');

% draw rois
r = cellArray(r);
overlayImage = nan(size(b.vtcs,1),3);
for rNum = 1:length(r)
  c = getSmoothColor(rNum,length(r),'hsv');
  if rNum ~= length(r)
    c = [1 1 0];
  else
    c = [1 0 0];
  end
  overlayImage(r{rNum}.vertices,1) = c(1);
  overlayImage(r{rNum}.vertices,2) = c(2);
  overlayImage(r{rNum}.vertices,3) = c(3);
end
patch('vertices', b.vtcs, 'faces', b.tris,'FaceVertexCData', overlayImage,'facecolor','interp','edgecolor','none','FaceAlpha',0.4);
axis equal;
axis off;
camva(4);
rotate3d on;
title(titleStr);
set(gca(f),'XDir','reverse');
set(gca(f),'YDir','normal');
set(gca(f),'ZDir','normal');
drawnow

%%%%%%%%%%%%%%%%%%%%%%%%
%    parseArguments    %
%%%%%%%%%%%%%%%%%%%%%%%%
function sessions = parseArguments(baseAnatomies,sessionPath,groupNum,scanNum,analysisName,overlayName,roiName)

sessions = [];

% check length of arguments, if we only have one, then copy for all scans
nSessions = length(sessionPath);
if nSessions <= 1
  disp(sprintf('(mlrComputeAtlasMean) Only %i sessions passed in',nSessions));
  return
end

% if an overlayName is passed in
if ~isempty(overlayName)
  % validate argument length of rest of arguments
  groupNum = validateArgumentLength('groupNum',groupNum,nSessions);
  if isempty(groupNum),return,end
  scanNum = validateArgumentLength('scanNum',scanNum,nSessions);
  if isempty(scanNum),return,end
  analysisName = validateArgumentLength('analysisName',analysisName,nSessions);
  if isempty(analysisName),return,end
  overlayName = validateArgumentLength('overlayName',overlayName,nSessions);
  if isempty(overlayName),return,end
  sessions.doOverlay = 1;
else
  sessions.doOverlay = 0;
end

% if an roi is passed in
if ~isempty(roiName)
  roiName = validateArgumentLength('roiName',roiName,nSessions);
  sessions.doROI = 1;
else
  sessions.doROI = 0;
end

% pack into output and return
sessions.nSessions = nSessions;
sessions.sessionPath = sessionPath;
if sessions.doOverlay
  sessions.groupNum = convertToCellArray(groupNum);
  sessions.scanNum = convertToCellArray(scanNum);
  sessions.analysisName = convertToCellArray(analysisName);
  sessions.overlayName = convertToCellArray(overlayName);
end
sessions.baseAnatomies = cellArray(baseAnatomies);
sessions.nBases = length(sessions.baseAnatomies);
if sessions.doROI
  sessions.roiName = cellArray(roiName);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    validateArgumentLength    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function arg = validateArgumentLength(argname,arg,nSessions)


% copy arg to all, if it is a single string that is passed in
if isstr(arg)
  thisArg = arg;
  arg = {};
  for i = 1:nSessions
    arg{i} = thisArg;
  end
elseif length(arg) == 1
  % arg is a number
  if isnumeric(arg)
    arg(2:nSessions) = arg(1);
  end
end

% check length
if length(arg) ~= nSessions
  disp(sprintf('(mlrComputeAtlasMean) %s length (%i) does not match number of sessions (%i)',argname,length(arg),nSessions));
  arg = [];
  return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    convertToCellArray    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a = convertToCellArray(a)

if ~iscell(a)
  for i = 1:length(a)
    b{i} = a(i);
  end
  a = b;
end

%%%%%%%%%%%%%%%%%
%    myisdir    %
%%%%%%%%%%%%%%%%%
function tf = myisdir(dirname)

tf = (length(dir(dirname)) ~= 0);
  