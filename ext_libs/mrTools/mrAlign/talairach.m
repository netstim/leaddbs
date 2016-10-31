% talairach.m
%
%        $Id$
%      usage: talinfo = talairach(volume)
%         by: justin gardner
%       date: 05/04/07
%    purpose: 
%             program to talairach a volume, call with filename
%             talinfo = talairach('jg041001');
%  
%             talinfo contains the talairach points. You can also
%             reset the taliarach points by doing:
%
%             talinfo = talairach(talinfo);
%

function talinfo = talairach(event,isEvent)

% check arguments
if ~any(nargin == [1 2 3])
    help talairach
    return
end

global gTalairach;

% init arguments
if nargin == 1
  % if we are passed in a structure then this should be a talinfo
  if isstruct(event)
    if ~isfield(event,'filename')
      disp(sprintf('(talairach) Input tal strucutre must have a filename'));
      return
    else
      talinfo = event;
      event = 'init';
    end
    % otherwise it is init event
  elseif isstr(event)
    talinfo.filename = event;
    event = 'init';
  end
end

switch (event)
  case 'init'
    talinfo = initHandler(talinfo);
  case 'end'
    endHandler;
  case 'mouseDown'
    mouseDownHandler;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mousedown
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mouseDownHandler

global gTalairach;

% get the pointer location
pointerLoc = get(gca,'CurrentPoint');
% get mouse, remembering that we have swapped y/x
% (see getImageSlice) 
mouseX = round(pointerLoc(1,1));
mouseY = round(pointerLoc(1,2));

% which figure we are on
if gcf == gTalairach.fig(1)
    a = [1 2 3];
    %  mouseX = size(gTalairach.vol,2)-mouseX+1;
    mouseY = size(gTalairach.vol,3)-mouseY+1;
    x = gTalairach.pos(1);
    y = mouseX;
    z = mouseY;
elseif gcf == gTalairach.fig(2)
    a = [2 1 3];
    %  mouseX = size(gTalairach.vol,1)-mouseX+1;
    mouseY = size(gTalairach.vol,3)-mouseY+1;
    x = mouseX;
    y = gTalairach.pos(2);
    z = mouseY;
elseif gcf == gTalairach.fig(3)
    a = [3 1 2];
    %  mouseX = size(gTalairach.vol,1)-mouseX+1;
    mouseY = size(gTalairach.vol,2)-mouseY+1;
    x = mouseX;
    y = mouseY;
    z = gTalairach.pos(3);
else
    return
end

% display the point being clicked on
%disp(sprintf('x:%i y:%i z:%i',x,y,z));

% set the position
if ((mouseX > 0) && (mouseX < gTalairach.dim(a(2))) && ...
    (mouseY > 0) && (mouseY < gTalairach.dim(a(3))))
    gTalairach.pos(a(2)) = mouseX;
    gTalairach.pos(a(3)) = mouseY;
    set(gTalairach.fig(a(1)),'pointer',mlrFullCrosshair);
else
    set(gTalairach.fig(a(1)),'pointer','arrow');
end

% refresh the display
refreshTalairachDisplay;

% and reset the current point int he controls
currentPos = xformView2Vol([x y z]);
params.currentX = currentPos(1);params.currentY = currentPos(2);params.currentZ = currentPos(3);
mrParamsSet(params);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   initTalairachControls %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initTalairachControls

global gTalairach;

paramsInfo = {};
paramsInfo{end+1} = {'setPoint',putOnTopOfList('None',gTalairach.refPoints),'Use this point for which talairach coordinate'};
paramsInfo{end+1} = {'currentX',gTalairach.pos(1),'incdec=[-1 1]',sprintf('minmax=[1 %i]',gTalairach.dim(1)),'The current selected point x coordinate'};
paramsInfo{end+1} = {'currentY',gTalairach.pos(2),'incdec=[-1 1]',sprintf('minmax=[1 %i]',gTalairach.dim(2)),'The current selected point y coordinate'};
paramsInfo{end+1} = {'currentZ',gTalairach.pos(3),'incdec=[-1 1]',sprintf('minmax=[1 %i]',gTalairach.dim(3)),'The current selected point z coordinate'};
for i = 1:length(gTalairach.refPoints)
  paramsInfo{end+1} = {gTalairach.refPoints{i},gTalairach.(gTalairach.refPoints{i}),sprintf('Talairach point %s',gTalairach.refPoints{i})};
end
% which plane to display
gTalairach.displayPlanes = {'None','AC-PC-SAC axis','AC-PC-SAC axial plane','AC-PC-SAC coronal plane','AC-PC-SAC sagittal plane'};
paramsInfo{end+1} = {'displayPlane',gTalairach.displayPlanes,'Which plane to display'};
paramsInfo{end+1} = {'displayPoints',1,'type=checkbox','Display point in yellow'};
paramsInfo{end+1} = {'straighten',{'None','Revert to original','AC-PC','AC-SAC','AC-PC-SAC'},'Rotate display to straighten across desired points. AC-PC will straighten in the sagittal and axial planes. AC-SAC will straighten only the coronal view. AC-PC-SAC applies both AC-PC and AC-SAC straightens.'};
paramsInfo{end+1} = {'showPoint',putOnTopOfList('None',gTalairach.refPoints),'Change to display to show selected point'};

% put up the dialog
params = mrParamsDialog(paramsInfo,'Set Talairach Points',[],@talairachControlsHandler,[],@endHandler,@cancelHandler);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   talairachControlsHandler %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function talairachControlsHandler(params)

global gTalairach;

currentPos = xformVol2View([params.currentX params.currentY params.currentZ]);
if any(currentPos ~= gTalairach.pos)
  gTalairach.pos(1) = currentPos(1);
  gTalairach.pos(2) = currentPos(2);
  gTalairach.pos(3) = currentPos(2);
end

% if the user selects to set this point
if ~strcmp(params.setPoint,'None')
  params.(params.setPoint)(1) = params.currentX;
  params.(params.setPoint)(2) = params.currentY;
  params.(params.setPoint)(3) = params.currentZ;
  params.setPoint = 'None';
  mrParamsSet(params);
end

% change the points in the global
for i = 1:length(gTalairach.refPoints)
  gTalairach.(gTalairach.refPoints{i}) = params.(gTalairach.refPoints{i});
end

% display a slice
if ~strcmp(params.displayPlane,'None')
  % set the coordinates of the ACPC slice
  gTalairach.ACPCSliceCoords = calcACPCSliceCoords(params.displayPlane);
else
  gTalairach.ACPCSliceCoords = [];
end

% set whether to display points
gTalairach.displayPoints = params.displayPoints;

% display asked for point
if ~strcmp(params.showPoint,'None')
  if ~any(gTalairach.(params.showPoint) <= 0)
    gTalairach.pos = xformVol2View(gTalairach.(params.showPoint));
    currentPos = gTalairach.(params.showPoint);
    params.currentX = currentPos(1);
    params.currentY = currentPos(2);
    params.currentZ = currentPos(3);
  end
  params.showPoint = 'None';
  mrParamsSet(params);
end

% straighten image if called for
if ~strcmp(params.straighten,'None')
gTalairach.vol2view = eye(4);
  % calculate the center of the volume xform
  centerOffset = [1 0 0 -gTalairach.dim(1)/2;...
		  0 1 0 -gTalairach.dim(2)/2;...
		  0 0 1 -gTalairach.dim(3)/2;...
		  0 0 0 1];

  % get points
  AC = xformVol2View(gTalairach.AC);
  PC = xformVol2View(gTalairach.PC);
  SAC = xformVol2View(gTalairach.SAC);
  
  % set the angle of rotation and make the appropriate rotation matrix
  if strcmp(params.straighten,'AC-PC')
    % get rotaion in sagittal plane
    rotMatrix = getRotMatrix(AC,PC,'xy')*getRotMatrix(AC,PC,'yz');
  elseif strcmp(params.straighten,'AC-SAC')
    rotMatrix = getRotMatrix(AC,SAC,'xz');
  elseif strcmp(params.straighten,'AC-PC-SAC')
    rotMatrix = getRotMatrix(AC,PC,'xy')*getRotMatrix(AC,PC,'yz')*getRotMatrix(AC,SAC,'xz');
  elseif strcmp(params.straighten,'Revert to original')
    rotMatrix = eye(4);
    gTalairach.vol2view = eye(4);
  end
  
  % get the xform matrix
  gTalairach.vol2view = inv(centerOffset) * inv(rotMatrix) * centerOffset;
  
  % xform volume
  gTalairach.viewVol = vol2viewVol(gTalairach.vol,inv(gTalairach.vol2view));

  % transform current position
  gTalairach.pos = xformVol2View([params.currentX params.currentY params.currentZ]);
  
  % set back to no straighten
  params.straighten = 'None';
  mrParamsSet(params);
end

refreshTalairachDisplay;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get current talairach rotation matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function viewVol = vol2viewVol(vol,view2vol);
% warp the volume
if ~isequal(view2vol,eye(4))
  disppercent(-inf,'Warping volume');
  swapXY = [0 1 0 0;1 0 0 0;0 0 1 0;0 0 0 1];
  % warpAffine3 uses yx, not xy
  viewVol = warpAffine3(vol,swapXY*view2vol*swapXY,nan,[],'linear');
  %viewVol = warpAffine3(vol,swapXY*view2vol*swapXY,nan,[],'nearest');
  disppercent(inf);
else
  viewVol = vol;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get current talairach rotation matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rotMatrix = getRotMatrix(point1,point2,plane)

% get the sides of the triangle for this rotation
switch (plane) 
  case {'xy'}
   opposite = point1(1) - point2(1);
   adjacent = point1(2) - point2(2);
  case {'yz'}
   opposite = point1(3) - point2(3);
   adjacent = point1(2) - point2(2);
  case {'xz'}
   opposite = point1(1) - point2(1);
   adjacent = point1(3) - point2(3);
end

% make sure we have a non-zero hypotenuse
hypotenuse = sqrt(opposite^2+adjacent^2);
if hypotenuse==0
  mrWarnDlg(sprintf('(talairach) %s angle of rotation invalid',plane));
  rotMatrix = eye(4);
  return
end

% and convert to a rotation matrix
switch (plane) 
  case {'xy'}
   angle = -asin(opposite/hypotenuse);
   c = cos(angle);
   s = sin(angle);
   rotMatrix = [c -s 0 0;s  c 0 0;0  0 1 0;0  0 0 1];
  case {'yz'}
  angle=asin(opposite/hypotenuse);
   c = cos(angle);
   s = sin(angle);
   rotMatrix = [1  0  0 0;0  c -s 0;0  s  c 0;0  0  0 1];
  case {'xz'}
   angle = asin(opposite/hypotenuse);
   c = cos(angle);
   s = sin(angle);
   rotMatrix = [c  0 -s 0;0  1  0 0;s  0  c 0;0  0  0  1];
end

% get cosine and sine of angle
disp(sprintf('(talairach) %s rotation = %0.1f deg',plane,angle*180/pi));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function endHandler(ok)

global gTalairach;

% only run this if it is not being called by someonw else
if ~gTalairach.shutdown
  gTalairach.shutdown = 1;
  for i = 1:3
    if ishandle(gTalairach.fig(i))
      close(gTalairach.fig(i));
    end
  end
  gTalairach.shutdown = 0;
  % set whether ok was pushed
  if nargin >= 1
    gTalairach.ok = ok;
  else
    gTalairach.ok = 1;
  end
  uiresume;
end

gTalairach.init = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cancelHandler

endHandler(0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init the interrogator handler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function talinfo = initHandler(talinfo)

global gTalairach;

% add img extension if it is not a single file nifti
if ~getext(talinfo.filename,'nii')
  talinfo.filename = sprintf('%s.img',stripext(talinfo.filename));
end
% check for file
if isfile(talinfo.filename)
  % read header
  hdr = mlrImageReadNiftiHeader(talinfo.filename);
  % check to see if this is an LPI volume
  if hdr.qform_code==1
    axisLabels = mlrImageGetAxisLabels(hdr.qform44);
    if ~strcmp(axisLabels.orient,'LPI')
      dispHeader
      mrWarnDlg(sprintf('(talairach) !!!! Volume %s is in %s orientation and should be LPI -> see mlrImageLoad to reorient volume !!!!',talinfo.filename,axisLabels.orient));
      dispHeader
    else
      disp(sprintf('(talairach) Volume %s is in LPI orientation',talinfo.filename));
    end
  end
  % if this is a 4D volume then we have to take a particular volume or
  % or the mean. Ask the user what to do.
  doMean = 0;subset = [];
  if hdr.dim(5) > 1
    paramsInfo = {{'volNum',0,'incdec=[-1 1]','round=1',sprintf('minmax=[0 %i]',hdr.dim(5)),'Choose volume number to display (0 for mean)'}};
    params = mrParamsDialog(paramsInfo,'Choose volume to use (0 for mean)');
    if isempty(params),return,end
    if params.volNum ~= 0
      subset = {[],[],[],params.volNum};
    else
      doMean = 1;
    end
  end
  % read it
  disppercent(-inf,sprintf('Loading %s',talinfo.filename));
  if ~isempty(subset)
    [vol hdr] = mlrImageReadNifti(talinfo.filename,subset);
  else
    [vol hdr] = mlrImageReadNifti(talinfo.filename);
  end
  if doMean,vol = mean(vol,4);end
  disppercent(inf);
else
  disp(sprintf('(talairach) Could not open file %s',talinfo.filename));
  return
end

% get figure handles
gTalairach = [];
gTalairach.fig(1) = mlrSmartfig('talairach1');
gTalairach.fig(2) = mlrSmartfig('talairach2');
gTalairach.fig(3) = mlrSmartfig('talairach3');
gTalairach.roiCoords = [];
gTalairach.shutdown = 0;
gTalairach.filename = talinfo.filename;
gTalairach.init = 1;
% set the callbacks appropriately
for i= 1:3
  % turn off menu bars 
  set(gTalairach.fig(i),'MenuBar','none');
  set(gTalairach.fig(i),'NumberTitle','off');
  % turn off close button
  set(get(gTalairach.fig(i),'Children'),'Visible','off')
  % set window button function
  set(gTalairach.fig(i),'WindowButtonDownFcn',sprintf('talairach(''mouseDown'',1)'));
  % and delete function
  set(gTalairach.fig(i),'DeleteFcn',sprintf('talairach(''end'',1)'));
end

% set pointer to crosshairs
set(gTalairach.fig(1),'pointer',mlrFullCrosshair);
set(gTalairach.fig(2),'pointer',mlrFullCrosshair);
set(gTalairach.fig(3),'pointer',mlrFullCrosshair);

% set volume
gTalairach.vol = vol;
gTalairach.hdr = hdr;
gTalairach.dim = size(vol);
gTalairach.gamma = 0.4;

gTalairach.pos(1) = round(size(vol,1)/2);
gTalairach.pos(2) = round(size(vol,2)/2);
gTalairach.pos(3) = round(size(vol,3)/2);

% compute color map
g = gray(256);
y = g;y(:,3) = 0;
c = g;c(:,1) = 0;
m = g;m(:,2) = 0;
r = g;r(:,2:3) = 0;
myColormap = [g;y;c;r;m;m];

gTalairach.params.width = 256;
gTalairach.params.height = 256;;
gTalairach.params.viewSlice = 1;
gTalairach.params.dispROI = 0;
gTalairach.params.xyRot = 0;
gTalairach.params.yzRot = 0;
gTalairach.params.xzRot = 0;
gTalairach.params.xCenter = 0;
gTalairach.params.yCenter = 0;
gTalairach.params.zCenter = 0;

% compute slice coordinates
gTalairach.ACPCSliceCoords = [];

% compute each dimension coordinates
[x y] = meshgrid(1:gTalairach.dim(2),1:gTalairach.dim(3));
gTalairach.sliceCoords{1}(1,1:length(x(:))) = gTalairach.pos(1);
gTalairach.sliceCoords{1}(2,:) = x(:);
gTalairach.sliceCoords{1}(3,:) = y(:);
gTalairach.sliceMin{1}(1:3) = [gTalairach.pos(1) min(x(:)) min(y(:))];
gTalairach.sliceMax{1}(1:3) = [gTalairach.pos(1) max(x(:)) max(y(:))];
[x y] = meshgrid(1:gTalairach.dim(1),1:gTalairach.dim(3));
gTalairach.sliceCoords{2}(1,:) = x(:);
gTalairach.sliceCoords{2}(2,:) = gTalairach.pos(2);
gTalairach.sliceCoords{2}(3,:) = y(:);
gTalairach.sliceMin{2}(1:3) = [min(x(:)) gTalairach.pos(1) min(y(:))];
gTalairach.sliceMax{2}(1:3) = [max(x(:)) gTalairach.pos(1) max(y(:))];
[x y] = meshgrid(1:gTalairach.dim(1),1:gTalairach.dim(2));
gTalairach.sliceCoords{3}(1,:) = x(:);
gTalairach.sliceCoords{3}(2,:) = y(:);
gTalairach.sliceCoords{3}(3,:) = gTalairach.pos(3);
gTalairach.sliceMin{3}(1:3) = [min(x(:)) min(y(:)) gTalairach.pos(1)];
gTalairach.sliceMax{3}(1:3) = [max(x(:)) max(y(:)) gTalairach.pos(1)];

% set up cache
gTalairach.c = mrCache('init',50);

% set the talairach ref points
gTalairach.refPoints = {'AC','PC','SAC','IAC','PPC','AAC','LAC','RAC'};
for i = 1:length(gTalairach.refPoints)
  if isfield(talinfo,gTalairach.refPoints{i})
    gTalairach.(gTalairach.refPoints{i}) = talinfo.(gTalairach.refPoints{i});
  else
    gTalairach.(gTalairach.refPoints{i}) = [0 0 0];
  end
end

% set to display the tal points
gTalairach.displayPoints = 1;

if ~isfield(talinfo,'vol2view')
  gTalairach.vol2view = eye(4);
else
  gTalairach.vol2view = talinfo.vol2view;
end

% transform the volume to the one that will be displayed,
% i.e. the "viewVol"
gTalairach.viewVol = vol2viewVol(gTalairach.vol,inv(gTalairach.vol2view));

% display three windows
figure(gTalairach.fig(1));
dispVolumeSlice(1,gTalairach.pos(1));
colormap(myColormap);

figure(gTalairach.fig(2));
dispVolumeSlice(2,gTalairach.pos(2));
colormap(myColormap);

figure(gTalairach.fig(3));
dispVolumeSlice(3,gTalairach.pos(3));
colormap(myColormap);

% set up controls handler
initTalairachControls;

% wait for ui to finish
uiwait;

% return talinfo properly
for i = 1:length(gTalairach.refPoints)
  talinfo.(gTalairach.refPoints{i}) = gTalairach.(gTalairach.refPoints{i});
end

% get the current view transformation matrix
talinfo.vol2view = gTalairach.vol2view;

% return empty if ok is not set
if ~gTalairach.ok
  talinfo = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% redisplay everything
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function refreshTalairachDisplay

global gTalairach;

% and redisplay
for i = 1:3
  figure(gTalairach.fig(i));
  dispVolumeSlice(i,gTalairach.pos(i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate ACPC slice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ACPCSliceCoords = calcACPCSliceCoords(slicePlane)

global gTalairach;

% calculate the AC shift xForm
ACOffset = [1 0 0 gTalairach.AC(1);...
	    0 1 0 gTalairach.AC(2);...
	    0 0 1 gTalairach.AC(3);...
	    0 0 0 1];

% get AC/PC
AC = gTalairach.AC;
PC = gTalairach.PC;
SAC = gTalairach.SAC;
dim = gTalairach.dim;

% get the initial coords
lineCoords = -dim(1):dim(1);
[x y] = meshgrid(lineCoords,lineCoords);
initCoords(1,:) = x(:);
initCoords(2,:) = y(:);
initCoords(3,:) = 0;
initCoords(4,:) = 1;

% compute coordinates after transformation
switch (slicePlane)
  case {'AC-PC-SAC axis'}
    zeroCoords = zeros(size(lineCoords));
    initCoords = [];
    initCoords(1,:) = [lineCoords zeroCoords zeroCoords];
    initCoords(2,:) = [zeroCoords lineCoords zeroCoords];
    initCoords(3,:) = [zeroCoords zeroCoords lineCoords];
    initCoords(4,:) = 1;
    ACPCSliceCoords = ACOffset*getRotMatrix(AC,PC,'xy')*getRotMatrix(AC,PC,'yz')*getRotMatrix(AC,SAC,'xz')*initCoords;
  case {'AC-PC-SAC axial plane'}
    ACPCSliceCoords = ACOffset*getRotMatrix(AC,PC,'xy')*getRotMatrix(AC,PC,'yz')*getRotMatrix(AC,SAC,'xz')*initCoords;
  case {'AC-PC-SAC coronal plane'}
    rotMatrix = [1 0 0 0;0 0 1 0;0 1 0 0;0 0 0 1];
    ACPCSliceCoords = ACOffset*getRotMatrix(AC,PC,'xy')*getRotMatrix(AC,PC,'yz')*getRotMatrix(AC,SAC,'xz')*rotMatrix*initCoords;
 case {'AC-PC-SAC sagittal plane'}
    rotMatrix = [0 0 1 0;0 1 0 0;1 0 0 0;0 0 0 1];
    ACPCSliceCoords = ACOffset*getRotMatrix(AC,PC,'xy')*getRotMatrix(AC,PC,'yz')*getRotMatrix(AC,SAC,'xz')*rotMatrix*initCoords;
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% display slice of volume
%%%%%%%%%%%%%%%%%%%%%%%%%
function slice = dispVolumeSlice(sliceDim,sliceNum)

global gTalairach

% save slice num
gTalairach.sliceNum(sliceDim) = sliceNum;

% get the other dimensions
otherDim = setdiff([1 2 3],sliceDim);

% clear axis so we don't keep drawing over old
cla

% get slice out of volume
switch (sliceDim)
  case {1}
    slice = squeeze(gTalairach.viewVol(sliceNum,:,:));
  case {2}
    slice = squeeze(gTalairach.viewVol(:,sliceNum,:));
  case {3}
    slice = squeeze(gTalairach.viewVol(:,:,sliceNum));
end

% update slice coords
gTalairach.sliceCoords{sliceDim}(sliceDim,:) = sliceNum;
gTalairach.sliceMin{sliceDim}(sliceDim) = sliceNum;
gTalairach.sliceMax{sliceDim}(sliceDim) = sliceNum;

% gamma correct
smax = max(slice(:));smin = min(slice(:));
slice = (slice-smin)./(smax-smin);
slice = 256*(slice.^gTalairach.gamma);

% now check to see if there is any overlap between this slice
% and the talairach
intersectCoords = [];displayIntersectCoords = [];
sMin = gTalairach.sliceMin{sliceDim};
sMax = gTalairach.sliceMax{sliceDim};

% this code used to use intersect, but that is hoeplessly
% slow. Because we are always displaying cardinal views
% we can just do min/max checking
if ~isempty(gTalairach.ACPCSliceCoords)
  ACPCSliceCoords = xformVol2View(gTalairach.ACPCSliceCoords(1:3,:)');
else
  ACPCSliceCoords = zeros(1,3);
end
thisSliceCoords = ...
    ((ACPCSliceCoords(:,1) >= sMin(1)) & ...
     (ACPCSliceCoords(:,1) <= sMax(1)) & ...
     (ACPCSliceCoords(:,2) >= sMin(2)) & ...
     (ACPCSliceCoords(:,2) <= sMax(2)) & ...
     (ACPCSliceCoords(:,3) >= sMin(3)) & ...
     (ACPCSliceCoords(:,3) <= sMax(3)));
displayIntersectCoords = ACPCSliceCoords(thisSliceCoords,:);
    
% get the indexes in the current slice for those intersection coordinates
if ~isempty(displayIntersectCoords)
    displayIntersectIndexes = sub2ind(gTalairach.dim(otherDim),displayIntersectCoords(:,otherDim(1)),displayIntersectCoords(:,otherDim(2)));
else
    displayIntersectIndexes = [];
end

% and those for the current slice to be cyan
slice(displayIntersectIndexes) = slice(displayIntersectIndexes)/2+640;

if isfield(gTalairach,'refPoints') && gTalairach.displayPoints
  for tnum = 1:length(gTalairach.refPoints)
    % now get the coords of each tal ref point
    talCoord = xformVol2View(gTalairach.(gTalairach.refPoints{tnum}));
    if ((talCoord(1) >= sMin(1)) & ...
	(talCoord(1) <= sMax(1)) & ...
	 (talCoord(2) >= sMin(2)) & ...
	(talCoord(2) <= sMax(2)) & ...
	(talCoord(3) >= sMin(3)) & ...
	(talCoord(3) <= sMax(3)))
      % get the indexes
      talIndexes = sub2ind(gTalairach.dim(otherDim),talCoord(otherDim(1)),talCoord(otherDim(2)));
      % and set them to yellow
      slice(talIndexes) = 512;
    end
  end
  % now get the current point
  currentPos = gTalairach.pos;
  if ((currentPos(1) >= sMin(1)) & ...
      (currentPos(1) <= sMax(1)) & ...
      (currentPos(2) >= sMin(2)) & ...
      (currentPos(2) <= sMax(2)) & ...
      (currentPos(3) >= sMin(3)) & ...
      (currentPos(3) <= sMax(3)))
    % get the indexes
    indexes = sub2ind(gTalairach.dim(otherDim),currentPos(otherDim(1)),currentPos(otherDim(2)));
    % and set it to red
    slice(indexes) = 1024;
  end
end

% flipud and transpose
slice = flipud(slice');

% display
image(slice);
axis off;axis square

titles = {sprintf('%s\nsagittal (left-right)',gTalairach.filename),sprintf('%s\ncoronal (back-front)',gTalairach.filename),sprintf('%s\naxial (bottom-top)',gTalairach.filename)};
title(titles{sliceDim},'Interpreter','none');
dimLabels = 'xyz';
xlabel(dimLabels(otherDim(1)));
ylabel(dimLabels(otherDim(2)));
% flip the y-axis up-down
set(gca,'YTickLabel',flipud(get(gca,'YTickLabel')));
YLim = get(gca,'YLim');
set(gca,'YTick',YLim(2)-fliplr(get(gca,'YTick')));
axis on
%%%%%%%%%%%%%%%%%%%%%%%%%%
% get talairach cache val
%%%%%%%%%%%%%%%%%%%%%%%%%%
function cacheID = getCacheID(sliceNum)

global gTalairach;
p = gTalairach.params;

cacheID = sprintf('%0.1f_%0.1f_%0.1f_%0.1f_%0.1f_%0.1f_%i_%i_%0.1f',p.yzRot,p.xzRot,p.xyRot,p.xCenter,p.yCenter,p.zCenter,p.width,p.height,sliceNum);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% xform from volume coordinates to the view coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%
function viewCoords = xformVol2View(volCoords)

global gTalairach;

if isequal(size(volCoords),[1 3])
  viewCoords = round(gTalairach.vol2view*[volCoords 1]');
else
  volCoords(:,4) = 1;
  viewCoords = round(gTalairach.vol2view*volCoords');
end
viewCoords = viewCoords(1:3,:)';

%%%%%%%%%%%%%%%%%%%%%%%%%%
% xform from volume coordinates to the view coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%
function volCoords = xformView2Vol(viewCoords)

global gTalairach;

if isequal(size(viewCoords),[1 3])
  volCoords = round(inv(gTalairach.vol2view)*[viewCoords 1]');
else
  volCoords = viewCoords;
end
