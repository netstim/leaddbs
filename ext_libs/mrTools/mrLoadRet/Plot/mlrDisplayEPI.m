% mlrDisplayEPI.m
%
%        $Id$
%      usage: mlrDisplayEPI(<v>,<groupNum>)
%         by: justin gardner
%       date: 08/16/07
%    purpose: display epis as a movie for inspection. 
%             Can be called with a view:
%             v = newView;
%             mlrDisplayEPI(v,'MotionComp')
%
%             or with a volume
%             data = mlrImageReadNifti('epiImages.hdr');
%             mlrDisplayEPI(data);
%
%             or with a filename
%             mlrDisplayEPI('epiImages.hdr');
%
function retval = mlrDisplayEPI(v,groupNum)

% check arguments
if ~any(nargin == [0 1 2])
  help mlrDisplayEPI
  return
end

global gMLRDisplayEPI;
gMLRDisplayEPI.deleteViewAtEnd = 0;
needsWarping = 0;
passedInVol = 0;

% setup a cache
gMLRDisplayEPI.c = mrCache('init',1000);
gMLRDisplayEPI.hdr = [];
      

% don't display too many frames at once
maxFramesAtOnce = 500;

% if we were passed in a strucutre
if ~ieNotDefined('v') && ~isview(v)
  % check for a file
  if isstr(v)
    filename = setext(v,'hdr');
    if isfile(filename)
      [v hdr] = mlrImageReadNifti(filename);
      if isempty(v),return,end
      gMLRDisplayEPI.hdr = hdr;
    else
      disp(sprintf('(mlrDisplayEPI) Could not open file %s',filename));
      return
    end
  end
  if isnumeric(v)
    % set up variables as appropriate
    gMLRDisplayEPI.v = v;
    gMLRDisplayEPI.nScans = 1;
    gMLRDisplayEPI.nSlices = size(v,3);
    gMLRDisplayEPI.nFrames = size(v,4);
    currentScan = 1;
    curSlice = 1;
    dialogTitle = sprintf('Displaying passed in volume');
    passedInVol = 1;
    clear v;
  else
    disp(sprintf('(mlrDisplayEPI) First argument must be a valid view, filename, or data'));
    return
  end
  gMLRDisplayEPI.fignum = mlrSmartfig('mlrDisplayEPI');
  gMLRDisplayEPI.graphWin = 0;
else
  % if no v, then create one
  if ieNotDefined('v'),
    % no view passed in, open a view
    v = newView;
    gMLRDisplayEPI.deleteViewAtEnd = 1;
  end
  % set the group
  if ~ieNotDefined('groupNum')
    v = viewSet(v,'curGroup',groupNum);
  end
  % keep v
  gMLRDisplayEPI.v = v;

  % get some parameters
  gMLRDisplayEPI.nScans = viewGet(v,'nScans');
  for iScan = 1:gMLRDisplayEPI.nScans
    gMLRDisplayEPI.nSlices(iScan) = viewGet(v,'nSlices',iScan);
    gMLRDisplayEPI.nFrames(iScan) = viewGet(v,'nFrames',iScan);
  end
  currentScan = viewGet(v,'currentScan');
  curSlice = viewGet(v,'curSlice');
  if isempty(curSlice),curSlice = 1;end
  dialogTitle = sprintf('EPI images for %s',viewGet(v,'groupName'));

  % get descriptions and scan2scan
  for scanNum = 1:viewGet(v,'nScans')
    descriptions{scanNum} = viewGet(v,'description',scanNum);
    scan2scan{scanNum} = viewGet(v,'scan2scan',1,[],scanNum);
    if ~isequal(scan2scan{scanNum},eye(4))
      needsWarping = 1;
    end
  end
  gMLRDisplayEPI.fignum = selectGraphWin([],'make new');
  gMLRDisplayEPI.graphWin = 1;
end

% set animating flag
gMLRDisplayEPI.animating = 0;
gMLRDisplayEPI.stopAnimating = 0;

% get interpMethod
gMLRDisplayEPI.interpMethod = mrGetPref('interpMethod');

% get max frames and slices
maxFrames = max(gMLRDisplayEPI.nFrames);
maxSlices = max(gMLRDisplayEPI.nSlices);

% set up params dialog
paramsInfo = {};
if ~passedInVol
  paramsInfo{end+1} = {'scanNum',currentScan,sprintf('minmax=[1 %i]',gMLRDisplayEPI.nScans),sprintf('incdec=[-1 1]'),'round=1','Choose scan to view','callback',@mlrDisplayEPIdeselectAll,'callbackArg=scan'};
  paramsInfo{end+1} = {'scanDescription',descriptions,'editable=0','group=scanNum','type=string','Description for scan'};
  paramsInfo{end+1} = {'scanNumMovie',0,'type=pushbutton','buttonString=Animate over scans','callback',@mlrDisplayEPIAnimate,'passParams=1','callbackArg','scanNum','Press to animate over scans'};
  paramsInfo{end+1} = {'displayAllScans',0,'type=checkbox','Click to display all scans at once','callback',@mlrDisplayEPIdeselectAll,'callbackArg=allExceptScan'};
else
  paramsInfo{end+1} = {'scanNum',1,'editable=0'};
end
paramsInfo{end+1} = {'sliceNum',curSlice,sprintf('minmax=[1 %i]',maxSlices),sprintf('incdec=[-1 1]'),'round=1','Choose slice number to view','callback',@mlrDisplayEPIdeselectAll,'callbackArg=slice'};
paramsInfo{end+1} = {'sliceNumMovie',0,'type=pushbutton','buttonString=Animate over slices','callback',@mlrDisplayEPIAnimate,'passParams=1','callbackArg','sliceNum','Press to animate over slices'};
paramsInfo{end+1} = {'displayAllSlices',0,'type=checkbox','Click to display all slices at once','callback',@mlrDisplayEPIdeselectAll,'callbackArg=allExceptSlice'};
paramsInfo{end+1} = {'frameNum',1,sprintf('minmax=[1 %i]',maxFrames),sprintf('incdec=[-1 1]'),'round=1','Choose frame to view','callback',@mlrDisplayEPIdeselectAll,'callbackArg=frame'};
paramsInfo{end+1} = {'frameNumMovie',0,'type=pushbutton','buttonString=Animate over frames','callback',@mlrDisplayEPIAnimate,'passParams=1','callbackArg','frameNum','Press to animate over frames'};
if all(gMLRDisplayEPI.nFrames < maxFramesAtOnce)
  paramsInfo{end+1} = {'displayAllFrames',0,'type=checkbox','Click to display all frames at once','callback',@mlrDisplayEPIdeselectAll,'callbackArg=allExceptFrame'};
end
if needsWarping
  paramsInfo{end+1} = {'warp',0,'type=checkbox','Apply warping implied by scan2scan transform to each image. This allows you to preview what will happen when you apply warping in averageTSeries or concatTSeries'};
  paramsInfo{end+1} = {'warpBaseScan',1,sprintf('minmax=[1 %i]',gMLRDisplayEPI.nScans),'incdec=[-1 1]','Choose which scan you want to warp the images to','contingent=warp'};
end
if ~isempty(gMLRDisplayEPI.hdr)
  paramsInfo{end+1} = {'displayNifitHdr',0,'type=pushbutton','buttonString=Display nifti header','callback',@mlrDisplayEPINiftiHeader};
end  

% display dialog
[gMLRDisplayEPI.f params] = mrParamsDialog(paramsInfo,dialogTitle,[],@mlrDisplayEPICallback,[],@mlrDisplayEPIClose);

% and draw first frame
mlrDisplayEPIDispImage(params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrDisplayEPINiftiHeader    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = mlrDisplayEPINiftiHeader(params)

global gMLRDisplayEPI;
hdr = gMLRDisplayEPI.hdr;
retval = 0;

numDims = hdr.dim(1);
dimstr = 'dim:';
pixdimstr = 'pixdim:';
for i = 1:numDims
  dimstr = sprintf('%s%i ',dimstr,hdr.dim(i+1));
  pixdimstr = sprintf('%s%i ',pixdimstr,hdr.pixdim(i+1));
end
disp(dimstr);
disp(pixdimstr);


% display qform and sform
disp(sprintf('++++++++++++++++++++++++++ qform ++++++++++++++++++++++++++'));
for rownum = 1:4
  disp(sprintf('%f\t%f\t%f\t%f',hdr.qform44(rownum,1),hdr.qform44(rownum,2),hdr.qform44(rownum,3),hdr.qform44(rownum,4)));
end
disp(sprintf('qformCode: %i',hdr.qform_code));

% display qform and sform
disp(sprintf('++++++++++++++++++++++++++ sform ++++++++++++++++++++++++++'));
for rownum = 1:4
  disp(sprintf('%f\t%f\t%f\t%f',hdr.sform44(rownum,1),hdr.sform44(rownum,2),hdr.sform44(rownum,3),hdr.sform44(rownum,4)));
end
disp(sprintf('sformCode: %i',hdr.sform_code));
  



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrDisplayEPIdeselectAll    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrDisplayEPIdeselectAll(deselectWhat,params)

if strcmp(deselectWhat,'slice')
  params.displayAllSlices = 0;
elseif strcmp(deselectWhat,'allExceptSlice')
  if isfield(params,'displayAllFrames'),params.displayAllFrames = 0;end
  if isfield(params,'displayAllScans'),params.displayAllScans = 0;end
elseif strcmp(deselectWhat,'frame')
  if isfield(params,'displayAllFrames'),params.displayAllFrames = 0;end
elseif strcmp(deselectWhat,'allExceptFrame')
  if isfield(params,'displayAllScans'),params.displayAllScans = 0;end
  params.displayAllSlices = 0;
elseif strcmp(deselectWhat,'scan')
  if isfield(params,'displayAllScans'),params.displayAllScans = 0;end
elseif strcmp(deselectWhat,'allExceptScan')
  if isfield(params,'displayAllFrames'),params.displayAllFrames = 0;end
  params.displayAllSlices = 0;
end
  
% turn off checkbox
mrParamsSet(params);
% and redisplay
mlrDisplayEPICallback(params);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  mlrDisplayEPIAnimate  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = mlrDisplayEPIAnimate(type,params)

retval = [];

global gMLRDisplayEPI;

% check to see whether we are already running an animation
if gMLRDisplayEPI.animating
  gMLRDisplayEPI.stopAnimating = 1;
  return
else
  gMLRDisplayEPI.animating = 1;
  gMLRDisplayEPI.stopAnimating = 0;
end

switch type
  case {'scanNum'}
   n = gMLRDisplayEPI.nScans;
   if isfield(params,'displayAllScans')
     params.displayAllScans = 0;
     mrParamsSet(params);
   end  
  case {'sliceNum'}
   params.displayAllSlices = 0;
   mrParamsSet(params);
   n = gMLRDisplayEPI.nSlices(params.scanNum);
  case {'frameNum'}
   if isfield(params,'displayAllFrames')
     params.displayAllFrames = 0;
     mrParamsSet(params);
   end
   n = gMLRDisplayEPI.nFrames(params.scanNum);
end

i = params.(type);
while (gMLRDisplayEPI.stopAnimating ~=1)
  % set the variable
  params.(type) = i;
  % set the dialog appropriately
  mrParamsSet(params);
  % load image
  mlrDisplayEPIDispImage(params);
  pause(.02);
  i = (mod(i,n))+1;
end

gMLRDisplayEPI.animating = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   mlrDisplayEPICallback   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrDisplayEPICallback(params)

global gMLRDisplayEPI;
if gMLRDisplayEPI.animating
  gMLRDisplayEPI.stopAnimating = 1;
end

% draw the image
mlrDisplayEPIDispImage(params);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   mlrDisplayEPIGetImage   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function epiImage = mlrDisplayEPIDispImage(params)

global gMLRDisplayEPI;
if ~isfield(params,'warp'),params.warp = 0;end
  
% check input arguments
nFrames = gMLRDisplayEPI.nFrames(params.scanNum);
nSlices = gMLRDisplayEPI.nSlices(params.scanNum);
params.frameNum = min(params.frameNum,nFrames);
params.sliceNum = min(params.sliceNum,nSlices);

% whether we have all slices or not
if params.displayAllSlices
  sliceNum = [1 nSlices];
else
  sliceNum = [params.sliceNum params.sliceNum];
end

% whether we have all frames or not
if isfield(params,'displayAllFrames') && params.displayAllFrames
  frameNum = [1 nFrames];
else
  params.displayAllFrames = 0;
  frameNum = [params.frameNum params.frameNum];
end

% check cache for image
if params.warp
  cacheStr = sprintf('%i_%s_%s_%i',params.scanNum,num2str(sliceNum),num2str(frameNum),params.warpBaseScan);
else
  cacheStr = sprintf('%i_%s_%s',params.scanNum,num2str(sliceNum),num2str(frameNum));
end
if isfield(params,'displayAllScans')
  cacheStr = sprintf('%s_%i',cacheStr,params.displayAllScans);
else
  params.displayAllScans = 0;
end

[epiImage gMLRDisplayEPI.c] = mrCache('find',gMLRDisplayEPI.c,cacheStr);

if isempty(epiImage)
  v = gMLRDisplayEPI.v;
  if isnumeric(v)
    epiImage = squeeze(v(:,:,[sliceNum(1):sliceNum(2)],[frameNum(1):frameNum(2)]));
  else
    % apply warping if necessary
    if params.warp
      % get scan2scan
      scan2scan = viewGet(v,'scan2scan',params.warpBaseScan,[],params.scanNum);

      if ~isequal(scan2scan,eye(4))
	% get dimensions of warp base scan - so we can transform
	% also if displaying all slices, we want to do this with
	% the number of slices in the warpBase
	warpBaseScanDims = viewGet(v,'scanDims',params.warpBaseScan);
	if params.displayAllSlices
	  sliceNum(2) = warpBaseScanDims(3);
	end
	
	% load the volume
	epiVolume = loadTSeries(v,params.scanNum,[],params.frameNum);

	% swapXY is needed because warpAffine3 is in yx not xy
	swapXY = [0 1 0 0;1 0 0 0;0 0 1 0; 0 0 0 1];

	% check to see if we also want to try to do a motion comp correction
	%if params.warpWithMotionComp
	%baseVolume = loadTSeries(v,params.warpBaseScan,[],params.frameNum);
	%scan2scan
	%  scan2scan = estMotionIter3(baseVolume,epiVolume,3,scan2scan)
	%end

	% compute transform
	M = swapXY * scan2scan * swapXY;

	% display transformation
	for rownum = 1:4
	  disp(sprintf('[%0.2f %0.2f %0.2f %0.2f]',M(rownum,1),M(rownum,2),M(rownum,3),M(rownum,4)));
	end
	disppercent(-inf,sprintf('Warping scan %i to match scan %i with transformation using %s',params.scanNum,params.warpBaseScan,gMLRDisplayEPI.interpMethod));
	epiVolume = warpAffine3(epiVolume,M,NaN,0,gMLRDisplayEPI.interpMethod,warpBaseScanDims);

	disppercent(inf);
	  epiImage = epiVolume(:,:,sliceNum(1):sliceNum(2));
      else
	% scan2scan was identity, so no warping is necessary
	epiImage = loadTSeries(v,params.scanNum,[],[frameNum(1) frameNum(2)]);
	epiImage = squeeze(epiImage(:,:,[sliceNum(1):sliceNum(2)],:));
      end
    else
      % no warping needed, just load from disk
      if params.displayAllScans
	% for alls cans, load each one
	epiImage = [];
	for scanNum = 1:gMLRDisplayEPI.nScans
	  if (gMLRDisplayEPI.nSlices(scanNum) >= sliceNum(1)) && (gMLRDisplayEPI.nFrames(scanNum) >= frameNum(1))
	    thisEpi = loadTSeries(v,scanNum,sliceNum(1),frameNum(1));
      sizeEpiImage = size(epiImage);
      if scanNum==1
        epiImage=thisEpi;
      elseif any(size(thisEpi)>sizeEpiImage(1:2)) %handle scans with different sizes
        epiImage2 = NaN([size(thisEpi) gMLRDisplayEPI.nScans]);
        epiImage2(1:sizeEpiImage(1),1:sizeEpiImage(2),1:scanNum-1) = epiImage;
        epiImage = epiImage2;
        epiImage(:,:,scanNum) = thisEpi;
      elseif any(size(thisEpi)<sizeEpiImage(1:2))
        sizeThisEpi = size(thisEpi);
        epiImage(:,:,scanNum)=NaN;
        epiImage(1:sizeThisEpi(1),1:sizeThisEpi(2),scnaNum) = thisEpi;
      else
        epiImage(:,:,scanNum) = thisEpi;
      end
	  end
	end
      else
	epiImage = loadTSeries(v,params.scanNum,[],[frameNum(1) frameNum(2)]);
	epiImage = squeeze(epiImage(:,:,[sliceNum(1):sliceNum(2)],:));
      end
    end
  end
    
  % Choose a sensible clipping value
  histThresh = length(epiImage(:))/1000;
  [cnt, val] = hist(epiImage(:),100);
  goodVals = find(cnt>histThresh);
  if isempty(goodVals)
    clipMin = 0;clipMax = 0;
  else
    clipMin = val(min(goodVals));
    clipMax = val(max(goodVals));
  end

  % and convert the image
  epiImage(epiImage<clipMin) = clipMin;
  epiImage(epiImage>clipMax) = clipMax;
  if (clipMax-clipMin) > 0
    epiImage = 255*(epiImage-clipMin)./(clipMax-clipMin);
  end
  % and save in cache
  gMLRDisplayEPI.c = mrCache('add',gMLRDisplayEPI.c,cacheStr,epiImage);
end

figure(gMLRDisplayEPI.fignum);

% now check if epiImage has multiple slices
nSlice = size(epiImage,3);
if nSlice > 1
  nRows = ceil(sqrt(nSlice));
  nCols = nRows;
  width = size(epiImage,2);
  height = size(epiImage,1);
  % make a giant image with all the slices
  newImage = zeros(height*nRows,width*nCols);
  k = 1;
  for i = 1:nRows
    for j = 1:nCols
      if (k <= nSlice)
	newImage((1+(i-1)*height):(1+i*height)-1,(1+(j-1)*width):(1+j*width)-1) = epiImage(:,:,k);
	k = k + 1;
      else
	break;
      end
    end
  end
  epiImage = newImage;
end

% always make it so that the h dimension is longer than v dimension
%if size(epiImage,1) > size(epiImage,2)
%  epiImage = epiImage';
%end

% display image
image(epiImage);
colormap(gray(256));
axis equal; axis tight; axis off
if params.displayAllSlices
  title(sprintf('Scan: %i slice 1:%i frame %i/%i',params.scanNum,nSlices,params.frameNum,nFrames));
else
  title(sprintf('Scan: %i slice %i/%i frame %i/%i',params.scanNum,params.sliceNum,nSlices,params.frameNum,nFrames));
end
drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%
%%   mlrDisplayEPIClose   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function mlrDisplayEPIClose

global gMLRDisplayEPI;

% delete the view if necessary
if gMLRDisplayEPI.deleteViewAtEnd
  deleteView(gMLRDisplayEPI.v);
end

% clear the global
gMLRDisplayEPI.stopAnimating = 1;
gMLRDisplayEPI.v = [];
gMLRDisplayEPI.c = [];

% close the graph window
if gMLRDisplayEPI.graphWin
  selectGraphWin([],'replace');
  closeGraphWin;
else
  close(gMLRDisplayEPI.fignum);
end

