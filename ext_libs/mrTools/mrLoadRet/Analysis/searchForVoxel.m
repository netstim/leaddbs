%
%      usage: searchForVoxel()
%         by: justin gardner
%       date: 11/08/07
%    purpose: interrogator function that displays desired voxel in
%             a different view
%
function retval = searchForVoxel(v,overlayNum,scan,x,y,s,roi,varargin)

% get the view number
viewNum = viewGet(v,'viewNum');

% get base types
for iBase = 1:viewGet(v,'numBase')
  baseType(iBase) = viewGet(v,'baseType',iBase);
end

% look for a non surface or flat map base
startBaseNum = first(find(baseType==0));
if ~isempty(startBaseNum)
  baseName = viewGet(v,'baseName',startBaseNum);
else
  baseName = viewGet(v,'baseName');
end

if ~isempty(viewGet(v,'currentOverlay')), hasOverlay = 1;else hasOverlay = 0;end

% get parameters
paramsInfo = {};
scanDims = viewGet(v,'scanDims');
sliceOrientations = {'default','axial','coronal','sagittal'};
baseNames = putOnTopOfList(baseName,viewGet(v,'baseNames'));
paramsInfo{end+1} = {'baseName',baseNames,'type=popupmenu','Choose which base to find voxel on'};
paramsInfo{end+1} = {'baseOrientation',sliceOrientations,'type=popupmenu','Choose which slice orientation to view voxel in'};
paramsInfo{end+1} = {'x',x,'incdec=[-1 1]','round=1',sprintf('minmax=[1 %i]',scanDims(1)),'Choose X coordinate to look for'};
paramsInfo{end+1} = {'y',y,'incdec=[-1 1]','round=1',sprintf('minmax=[1 %i]',scanDims(2)),'Choose Y coordinate to look for'};
paramsInfo{end+1} = {'s',s,'incdec=[-1 1]','round=1',sprintf('minmax=[1 %i]',scanDims(3)),'Choose S coordinate to look for'};
paramsInfo{end+1} = {'color',color2RGB,'type=popupmenu','Choose color to display voxel in'};
paramsInfo{end+1} = {'maxSearchRadius',5,'type=numeric','If there is no exact match, then will display the closest voxel within this search radius. Set smaller to force only display of more exact matches. Set higher if you want to allow matches that are farther away.'};
paramsInfo{end+1} = {'continuousMode',~isempty(startBaseNum),'type=checkbox','Opens another window with your anatomy that will be continuously update with the position as you move your mouse over points in the viewer','callback',@continuousModeCallback,'passParams=1','callbackArg',v,'editable',isempty(startBaseNum)};
if hasOverlay
  paramsInfo{end+1} = {'dispOverlay',0,'type=checkbox','Sets whether to display the current overlay in continuousMode','contingent=continuousMode'};
  paramsInfo{end+1} = {'overlayAlpha',1,'type=numeric','minmax=[0 1]','incdec=[-0.1 0.1]','Sets whether to display the current overlay in continuousMode','contingent=continuousMode'};
end

% start continuous mode
if ~isempty(startBaseNum)
  continuousModeCallback(v,mrParamsDefault(paramsInfo));
end

% put up dialog box
params = mrParamsDialog(paramsInfo,'Look for voxel','callback',@searchForVoxelCallback,'modal=1');

% shut down the continuousModeMouseMove function
continuousModeEnd;

if isempty(params)
  % if user hit cancel, then refresh
  % display to remove any marked voxel and return
  refreshMLRDisplay(viewNum);
  return
end

% switch to the chosen base
if ~strcmp(viewGet(v,'baseName'),params.baseName)
  v = viewSet(v,'curBase',viewGet(v,'baseNum',params.baseName));
end

% now switch to the appropriate orientation and slice
baseType = viewGet(v,'baseType');
if baseType == 0
  % check orientation
  whichSliceOrientation = find(strcmp(params.baseOrientation,sliceOrientations));
  % if user didn't ask for the default, set orientation apporpriately
  if whichSliceOrientation ~= 1
    v = viewSet(v,'sliceOrientation',whichSliceOrientation-1);
  end
  % first get which dimension we are looking at
  baseSliceIndex = viewGet(v,'baseSliceIndex');
  % now find voxel in base coordinates
  scan2base = inv(viewGet(v,'base2scan'));
  baseVoxel = round(scan2base*[params.x params.y params.s 1]');
  % now get the slice to switch to
  sliceNum = baseVoxel(baseSliceIndex);
  % make sure the asked for slice exists
  baseDims = viewGet(v,'baseDims');
  if ((sliceNum > 0) && (sliceNum <= baseDims(baseSliceIndex)))
    % if so, switch to it
    viewSet(v,'curSlice',sliceNum);
  end
end

% refresh the display
refreshMLRDisplay(viewNum);
hold on

% get the refreshed view
v = viewGet([],'view',viewNum);
  
% and then display the voxel
if baseType <= 1
  % find voxel in base coordinates
  scan2base = inv(viewGet(v,'base2scan'));
  baseVoxel = round(scan2base*[params.x params.y params.s 1]');

  % find the coordinates of the view
  baseCoords = round(viewGet(v,'cursliceBaseCoords'));

  % compute distance to matching voxel
  matchDistance = sqrt((baseCoords(:,:,1) - baseVoxel(1)).^2 + (baseCoords(:,:,2) - baseVoxel(2)).^2 + (baseCoords(:,:,3) - baseVoxel(3)).^2);
  minMatchDistance = min(matchDistance(:));
  
  % if there is no match within 5 voxels, give up
  if minMatchDistance > params.maxSearchRadius
    mrWarnDlg(sprintf('(searchForVoxel) No voxel within %f of searched for voxel. Closest voxel is within a radius of %f',params.maxSearchRadius,minMatchDistance));
    hold off
    return
  end

  % display info
  disp(sprintf('(searchForVoxel) Scan voxel [%i %i %i] projects to [%i %i %i] in %s. Displaying %i voxels within a radius of %f.',params.x,params.y,params.s,baseVoxel(1),baseVoxel(2),baseVoxel(3),params.baseName,length(find(matchDistance==minMatchDistance)),minMatchDistance));

  [y x] = find(matchDistance==minMatchDistance);
  for i = 1:length(x)
    plot(x,y,'.','Color',color2RGB(params.color));
  end
else
  % if we need to display on a surface,
  % first get base coord of point
  scan2base = inv(viewGet(v,'base2scan'));
  baseVoxel = round(scan2base*[params.x params.y params.s 1]');
  
  % and plot
  plot3(baseVoxel(1),baseVoxel(2),baseVoxel(3),'.','Color',color2RGB(params.color),'MarkerSize',30);
end

hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    continuousModeCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = continuousModeCallback(v,params)

retval = [];

if params.continuousMode
  continuousModeStart(v,params);
else
  continuousModeEnd;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   continuousModeEnd   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function continuousModeEnd(x,y,z)

global gCMode;

if isfield(gCMode,'continuousMode') && gCMode.continuousMode
  v = viewGet([],'view',gCMode.viewNum);
  if isempty(v),return,end
  fignum = viewGet(v,'fignum');
  set(fignum,'WindowButtonMotionFcn',gCMode.oldMouseMoveFcn);
  set(fignum,'WindowButtonDownFcn',gCMode.oldMouseDownFcn);
  close(gCMode.graphFig);
  gCMode.continuousMode = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   continuousModeStart   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function continuousModeStart(v,params)

% get some info that will be used by continuousModeMouseMove
global gCMode;
if ~isfield(gCMode,'continuousMode') || ~gCMode.continuousMode
  gCMode.fignum = viewGet(v,'fignum');
  g = guidata(gCMode.fignum);
  gCMode.axesNum = g.axis;
  gCMode.viewNum = viewGet(v,'viewNum');;
  
  % set the parameters
  updateContinuousModeParams(params);

  % set the main viewer to use our mouse move function
  gCMode.oldMouseMoveFcn = get(gCMode.fignum,'WindowButtonMotionFcn');
  gCMode.oldMouseDownFcn = get(gCMode.fignum,'WindowButtonDownFcn');
  set(gCMode.fignum,'WindowButtonMotionFcn',@continuousModeMouseMove);
  set(gCMode.fignum,'WindowButtonDownFcn',@continuousModeEnd);
  
  % open the figure for displaying
  gCMode.graphFig = selectGraphWin;
  subplot(1,3,1)
  gCMode.a1 = gca;
  axis(gCMode.a1,'off');
  subplot(1,3,2)
  gCMode.a2 = gca;
  axis(gCMode.a2,'off');
  subplot(1,3,3)
  gCMode.a3 = gca;
  axis(gCMode.a3,'off');
  colormap(gray);
  gCMode.continuousMode = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   updateContinuousModeParams   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updateContinuousModeParams(params)

global gCMode;

if isfield(gCMode,'viewNum')
  v = viewGet([],'view',gCMode.viewNum);
  gCMode.baseName = params.baseName;
  gCMode.baseNum = viewGet(v,'baseNum',params.baseName);
  gCMode.base2base = viewGet(v,'base2base',gCMode.baseNum);
  gCMode.baseType = viewGet(v,'baseType',gCMode.baseNum);
  gCMode.color = color2RGB(params.color);
  gCMode.baseOrientation = find(strcmp(params.baseOrientation,{'default','axial','coronal','sagittal'}));
  gCMode.base = viewGet(v,'baseData',gCMode.baseNum);
  gCMode.baseDims = size(gCMode.base);
  gCMode.base2overlay = viewGet(v,'base2scan',[],[],gCMode.baseNum);
  if isfield(params,'dispOverlay')
    gCMode.dispOverlay = params.dispOverlay;
    gCMode.overlayAlpha = params.overlayAlpha;
  else
    gCMode.dispOverlay = 0;
    gCMode.overlayAlpha = 0;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    continousModeMouseMove    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function continuousModeMouseMove(fignum,x);

global gCMode;
if ~isfield(gCMode,'continuousMode') || (gCMode.continuousMode==0)
  return
end
if (gCMode.baseType ~= 0) || isempty(gCMode.base2base)
  cla(gCMode.a1);
  axis(gCMode.a1,'off');
  if isempty(gCMode.base2base)
    title(sprintf('Could not compute xform to selected base - has base been aligned?'),'Parent',gCMode.a1);
  else
    title(sprintf('Continuous mode not implemented for surface/flat maps'),'Parent',gCMode.a1);
  end
  cla(gCMode.a2);
  axis(gCMode.a2,'off');
  cla(gCMode.a3);
  axis(gCMode.a3,'off');
  return
end
% get the mouse coordinates
[x y s xBase yBase sBase] = getMouseCoords;

% get the view
v = viewGet([],'view',gCMode.viewNum);

% convert to the base units
if ~isnan(xBase)
  baseCoord = gCMode.base2base * [xBase yBase sBase 1]';
  if (baseCoord(1) >= 1) && (baseCoord(1) <= gCMode.baseDims(1)) && (baseCoord(2) >= 1) && (baseCoord(2) <= gCMode.baseDims(2)) && (baseCoord(3) >= 1) && (baseCoord(3) <= gCMode.baseDims(3))
    if gCMode.baseOrientation == 1
      % clear the axis
      cla(gCMode.a1,'reset');
      % display the fitrst image
      if gCMode.dispOverlay
	% get base image
	[base.im,base.coords,base.coordsHomogeneous] = getBaseSlice(v,round(baseCoord(3)),3,0,gCMode.baseNum,gCMode.baseType);
	baseMin = min(base.im(:));baseMax = max(base.im(:));
	baseRGB = repmat((base.im-baseMin)./(baseMax-baseMin),[1 1 3]);
	% get overlay
	o = computeOverlay(v,gCMode.base2overlay,base.coordsHomogeneous,[gCMode.baseDims(1:2) 1],1);
	% composite
	img = (1-gCMode.overlayAlpha*o.alphaMap).*baseRGB + gCMode.overlayAlpha*o.alphaMap.*o.RGB;
	% and display
	image(img,'Parent',gCMode.a1);
      else
	imagesc(squeeze(gCMode.base(:,:,round(baseCoord(3)))),'Parent',gCMode.a1);
      end
      % display the current position
      hold(gCMode.a1,'on');
      axis(gCMode.a1,'equal');
      plot(baseCoord(2),baseCoord(1),'.','MarkerEdgeColor',gCMode.color,'Parent',gCMode.a1);
      title(sprintf('BaseCoord: %i %i %i',round(baseCoord(1)),round(baseCoord(2)),round(baseCoord(3))),'Parent',gCMode.a1);
      % display the second image
      cla(gCMode.a2,'reset');
      if gCMode.dispOverlay
	[base.im,base.coords,base.coordsHomogeneous] = getBaseSlice(v,round(baseCoord(2)),2,0,gCMode.baseNum,gCMode.baseType);
	baseMin = min(base.im(:));baseMax = max(base.im(:));
	baseRGB = repmat((base.im-baseMin)./(baseMax-baseMin),[1 1 3]);
	o = computeOverlay(v,gCMode.base2overlay,base.coordsHomogeneous,[gCMode.baseDims([1 3]) 1],1);
	img = (1-gCMode.overlayAlpha*o.alphaMap).*baseRGB + gCMode.overlayAlpha*o.alphaMap.*o.RGB;
	image(img,'Parent',gCMode.a2);
      else
	imagesc(squeeze(gCMode.base(:,round(baseCoord(2)),:)),'Parent',gCMode.a2);
      end
      % display the current position
      hold(gCMode.a2,'on');
      axis(gCMode.a2,'equal');
      plot(baseCoord(3),baseCoord(1),'.','MarkerEdgeColor',gCMode.color,'Parent',gCMode.a2);


      %display the third image
      cla(gCMode.a3,'reset');
      if gCMode.dispOverlay
	[base.im,base.coords,base.coordsHomogeneous] = getBaseSlice(v,round(baseCoord(1)),1,0,gCMode.baseNum,gCMode.baseType);
	baseMin = min(base.im(:));baseMax = max(base.im(:));
	baseRGB = repmat((base.im-baseMin)./(baseMax-baseMin),[1 1 3]);
	o = computeOverlay(v,gCMode.base2overlay,base.coordsHomogeneous,[gCMode.baseDims(2:3) 1],1);
	img = (1-gCMode.overlayAlpha*o.alphaMap).*baseRGB + gCMode.overlayAlpha*o.alphaMap.*o.RGB;
	image(img,'Parent',gCMode.a3);
      else
	imagesc(squeeze(gCMode.base(round(baseCoord(1)),:,:)),'Parent',gCMode.a3);
      end
      % display the current position
      hold(gCMode.a3,'on');
      axis(gCMode.a3,'equal');
      plot(baseCoord(3),baseCoord(2),'.','MarkerEdgeColor',gCMode.color,'Parent',gCMode.a3);
    end
  else
    axis(gCMode.a1,'off');
    title(sprintf('Voxel: [%s] not found in base %s with dims [%s]',num2str([xBase yBase sBase],'% 0.0f'),gCMode.baseName,num2str(gCMode.baseDims,' %i')),'Parent',gCMode.a1);
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get current mouse position in image coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x y s xBase yBase sBase xTal yTal zTal] = getMouseCoords;

global gCMode;
fignum = gCMode.fignum;
axesNum = gCMode.axesNum;
v = viewGet([],'view',gCMode.viewNum);

% no bases
if viewGet(v,'numBase') == 0
  x = nan;y= nan;s = nan;xBase = nan;yBase = nan;sBase = nan;xTal = nan;yTal = nan;zTal = nan;
  return
end
baseType = viewGet(v,'baseType');

% get location of pointer
pointerLoc = get(axesNum,'CurrentPoint');

if baseType <= 1
  mouseY = round(pointerLoc(1,1));
  mouseX = round(pointerLoc(1,2));

  % get base coordinates
  baseCoords = viewGet(v,'cursliceBaseCoords');
  % convert mouse to baseCoords
  if (mouseX>0) && (mouseX<=size(baseCoords,1)) && (mouseY>0) && (mouseY<=size(baseCoords,2))
    xBase = baseCoords(mouseX,mouseY,1);
    yBase = baseCoords(mouseX,mouseY,2);
    sBase = baseCoords(mouseX,mouseY,3);
  else
    x = nan;y = nan; s = nan;
    xBase = nan;yBase = nan; sBase = nan;
    xTal = nan; yTal = nan; zTal = nan;
    return
  end
else
  % handle getting coordinates for surface
  baseSurface = viewGet(v,'baseSurface');
  baseDims = viewGet(v,'baseSurfaceDims');
  pos = [];xBase = nan; yBase = nan; sBase = nan;
  % check mouse bounding box coords against baseDims
  % for a quick check to see if we are in the volume
  if all(pointerLoc(1,:) >= 0)
    % then use select3d which is slooow, but accurate
    hobj = get(axesNum,'Children');
    % make sure we are using the correct object (should be the 3D
    % brain). Use end here because with searchForVoxel we plot a
    % point on the image and the brain object always seems to be
    % last. But if this is not always the case, then we may need
    % to do a little more work here to find the correct object
    [pos vertex vi] = select3d(hobj(end));
    % convert the index to the coordinates
    if ~isempty(pos)
      baseCoordMap = viewGet(v,'baseCoordMap');
      %we'll take the coordinates of the middle of whatever range of cortical depth is currenlty selected
      corticalSlice = ceil(mean(viewGet(v,'corticalDepth'))*size(baseCoordMap.coords,5));
      pos = round(squeeze(baseCoordMap.coords(1,vi,1,:,corticalSlice)));
      xBase = pos(1);yBase = pos(2);sBase = pos(3);
    end
  end      
end

% transform from base coordinates to talairach coordinates, if the base has a talairach transform defined
base2tal = viewGet(v,'base2tal'); % keyboard
if(~isempty(base2tal))
  talCoords = round(base2tal * [xBase yBase sBase 1]');
  xTal = talCoords(1); yTal = talCoords(2); zTal = talCoords(3);
else
  xTal = nan; yTal = nan; zTal = nan;
end


% transform from base coordinates into scan coordinates
base2scan = viewGet(v,'base2scan');
if isempty(base2scan), x = nan; y=nan; s=nan; return,  end
transformed = round(base2scan*[xBase yBase sBase 1]');

x = transformed(1);
y = transformed(2);
s = transformed(3);

% get the scan dims to make sure we haven't jumped off end
% of scan
scanDims = viewGet(v,'scanDims',viewGet(v,'curScan'));
if ((x < 1) || (x > scanDims(1)) || ...
        (y < 1) || (y > scanDims(2)) || ...
        (s < 1) || (s > scanDims(3)))
    x = nan;y = nan;s = nan;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   searchForVoxelCallback   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function params = searchForVoxelCallback(params)

updateContinuousModeParams(params)
