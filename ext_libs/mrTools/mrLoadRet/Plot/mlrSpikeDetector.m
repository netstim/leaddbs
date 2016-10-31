% spikedetector.m
%
%        $Id$	
%      usage: spikeInfo = mlrSpikeDetector(scanNum,groupNum,<'criterion=10'>,<'dispfigs=1'>
%         by: justin gardner
%       date: 01/09/06
%    purpose: looks for spiking artifact in data set
%       e.g.: spikedetector(1,3)
% 
%             criterion sets how many std above mean a fourier
%             component has to be to be considered a
%             spike. (default=10)
%
function v = mlrSpikeDetector(v,scanNum,groupNum,varargin)

% check arguments
if nargin < 3
  help mlrSpikeDetector
  return
end

% get optional arguments
eval(evalargs(varargin));

if exist('criterion')~=1,criterion = 4;,end
if exist('dispfigs')~=1,dispfigs = 1;,end
if exist('recompute')~=1,recompute = 0;,end

% try to load spikeInfo from view
spikeInfo = viewGet(v,'spikeInfo',scanNum,groupNum);

if isempty(spikeInfo) || recompute
  % ask user how to recalculate
  paramsInfo{1} = {'criterion',criterion,'incdec=[-1 1]','minmax=[0 inf]','Criterion for spike detection. Spike detection works by computing the mean and standard deviation of each fourier component of the images across time. If a fourier component on any single volume exceeds criterion standard deviations of the mean, it is considered to be a spike. Default value is 4. i.e. A fourier component has to be 4 standard deviations greater from the mean to be considered to be a spike. This is a low threshold, but can be changed later on'};
  paramsInfo{2} = {'useMedian', 0, 'type=checkbox', 'Use the median and interquartile range to calculate the center and spread of the data.  This is useful if the data are very noisy.'};
  paramsInfo = mrParamsDialogSelectScans(v,groupNum,paramsInfo,scanNum);
  params = mrParamsDialog(paramsInfo,'mlrSpikeDetector params');
  if isempty(params)||~any(params.include),return,end
  % go through and recalculate
  scanNums = find(params.include);
  for s = scanNums
    spikeInfo = calcSpikeInfo(v,s,groupNum,params);
    v = viewSet(v,'spikeInfo',spikeInfo,s,groupNum);
  end
  saveSession;
  % now get back the one we want
  if ~ismember(scanNum,scanNums)
    scanNum = scanNums(1);
  end
  spikeInfo = viewGet(v,'spikeInfo',scanNum,groupNum);
end
if dispfigs
  hFigure = selectGraphWin;
  gData.spikeInfo=spikeInfo;
  gData.v = v;
  guidata(hFigure,gData);
  initSpikeFigure(hFigure);
  spikePlotController(hFigure);
end

%%%%%%%%%%%%%%%%%%%%%%%
%%   calcSpikeInfo   %%
%%%%%%%%%%%%%%%%%%%%%%%
function spikeInfo = calcSpikeInfo(v,scanNum,groupNum,params)

% load TSeries
spikeInfo.scanNum = scanNum;
spikeInfo.groupNum = groupNum;
spikeInfo.filename = viewGet(v,'tSeriesFile',scanNum,groupNum);
if isempty(spikeInfo.filename)
  disp(sprintf('(mlrSpikeDetector) Could not find scan %i in group %i',scanNum,groupNum));
  spikeInfo = [];
  return
end
  
% load file
v = viewSet(v,'curGroup',groupNum);
disppercent(-inf,sprintf('(mlrSpikeDetector) Loading time series for scan %i, group %i',scanNum,groupNum));
data = loadTSeries(v,scanNum);
% Dump junk frames
junkFrames = viewGet(v, 'junkframes', scanNum);
nFrames = viewGet(v, 'nFrames', scanNum);
data = data(:,:,:,junkFrames+1:junkFrames+nFrames);

spikeInfo.dim = size(data);
disppercent(inf);

% compute timecourse means for later 
for slicenum = 1:spikeInfo.dim(3)
  % calculate mean timecourse for slice
  sliceMeans(slicenum,:) = nanmean(nanmean(data(:,:,slicenum,:),1),2);
  % normalize to % signal change
  sliceMeans(slicenum,:) = 100*sliceMeans(slicenum,:)/mean(sliceMeans(slicenum,:));
end

% compute fourier transform of data
% calculating fourier transform of data
disppercent(-inf,'(mlrSpikeDetector) Calculating FFT');
% skip some frames in the beginning to account
% for saturation
if junkFrames < 5
    startframe = min(5,spikeInfo.dim(4));
else
    startframe = 1;
end
data = data(:,:,:,startframe:spikeInfo.dim(4));
for i = 1:size(data,4)
  disppercent(i/size(data,4));
  for j = 1:size(data,3)
    %first need to remove NaNs from data
    %let's replace them by the mean of each image
    thisData =  data(:,:,j,i);
    thisData(isnan(thisData)) = nanmean(thisData(:));
    %compute the spatial fourier transform
    data(:,:,j,i) = abs(fftshift(fft2(thisData)));
  end
end
disppercent(inf);

% get mean and std
if params.useMedian
    disppercent(-inf,'(mlrSpikeDetector) Calculating median and iqr');
    for slicenum = 1:spikeInfo.dim(3)
        disppercent(slicenum/spikeInfo.dim(3));
        meandata(:,:,slicenum) = squeeze(median(data(:,:,slicenum,:),4));
        stddata(:,:,slicenum) = squeeze(iqr(data(:,:,slicenum,:),4));
    end
else
    disppercent(-inf,'(mlrSpikeDetector) Calculating mean and std');
    for slicenum = 1:spikeInfo.dim(3)
        meandata(:,:,slicenum) = squeeze(mean(data(:,:,slicenum,:),4));
        stddata(:,:,slicenum) = squeeze(std(data(:,:,slicenum,:),0,4));
    end
end
disppercent(inf);


% now subtract off mean and see
% if there are any points above std criterion
slice = [];time = [];numspikes = [];spikelocs = {};meanZvalue=[];
disppercent(-inf,'(mlrSpikeDetector) Looking for spikes');
for i = 1:size(data,4)
  disppercent(i/spikeInfo.dim(4));
  data(:,:,:,i) = squeeze(data(:,:,:,i))-meandata;
  % see if any voxels are larger then expected
  for slicenum = 1:spikeInfo.dim(3)
    [spikex spikey] = find(squeeze(data(:,:,slicenum,i)) > params.criterion*squeeze(stddata(:,:,slicenum)));
    if ~isempty(spikex)
      slice(end+1) = slicenum;
      time(end+1) = startframe + i -1;
      numspikes(end+1) = length(spikex);
      spikelocs{end+1}.x = spikex;
      spikelocs{end}.y = spikey;
      spikelocs{end}.linear = find(squeeze(data(:,:,slicenum,i)) > params.criterion*squeeze(stddata(:,:,slicenum)));
      for iSpike = 1:length(spikex)
        spikelocs{end}.zValue(iSpike) =  data(spikex(iSpike),spikey(iSpike),slicenum,i)/squeeze(stddata(spikex(iSpike),spikey(iSpike),slicenum));
      end
      meanZvalue(end+1) = mean(spikelocs{end}.zValue);
    end
  end
end
disppercent(inf);
if length(slice)
  disp(sprintf('======================================================'));
  disp(sprintf('(mlrSpikeDetector) Found %i spikes at z>%.2f in scan %i, group %i',length(slice),params.criterion,scanNum,groupNum));
  disp(sprintf('======================================================'));
else
  disp(sprintf('(mlrSpikeDetector) No spikes in scan %i, group %i',scanNum,groupNum));
end

% pass them back in spikeInfo
spikeInfo.n = length(slice);
spikeInfo.slice = slice;
spikeInfo.time = time;
spikeInfo.numspikes = numspikes;
spikeInfo.spikelocs = spikelocs;
spikeInfo.criterion = params.criterion;
spikeInfo.sliceMeans = single(sliceMeans);
spikeInfo.meanZvalue = meanZvalue;
spikeInfo.maxZvalue = [zeros(spikeInfo.dim(3),startframe-1) permute(max(max(data./repmat(stddata,[1 1 1 size(data,4)]),[],1),[],2),[3 4 1 2])];

%%%%%%%%%%%%%%%%%%%%%%%%%
%%   initSpikeFigure   %%
%%%%%%%%%%%%%%%%%%%%%%%%%
function initSpikeFigure(hFigure)

gData = guidata(hFigure);

if isempty(gData.spikeInfo)
  return
end
  
% plot the slice means
gData.colorOrder = hsv(gData.spikeInfo.dim(3)+2);
gData.colorOrder = gData.colorOrder(3:end,:);
set(hFigure,'name',...
  sprintf('Spike Detector - %s:%i %s (%s)',...
           viewGet(gData.v,'groupName',gData.spikeInfo.groupNum),gData.spikeInfo.scanNum,...
           viewGet(gData.v,'description',gData.spikeInfo.scanNum,gData.spikeInfo.groupNum),...
           gData.spikeInfo.filename));

gData.nrows = 4;
gData.ncols = 5;
plot1 = 1:3;
plot2 = [6:8,11:13,16:18];
gData.hImagePlot = subplot(gData.nrows,gData.ncols,[4 5 9 10],'parent',hFigure);
gData.hFftPlot = subplot(gData.nrows,gData.ncols,[14 15 19 20],'parent',hFigure);

gData.hTseriesAxis = subplot(gData.nrows,gData.ncols,plot1,'parent',hFigure);
hold(gData.hTseriesAxis,'on');
set(gData.hTseriesAxis,'xLim',[0 gData.spikeInfo.dim(4)]);
gData.tSeriesYlim = [floor(min(gData.spikeInfo.sliceMeans(:))) ceil(max(gData.spikeInfo.sliceMeans(:)))];
set(gData.hTseriesAxis,'yLim',gData.tSeriesYlim);
gData.hCursorTseries = plot(gData.hTseriesAxis,[1 1],gData.tSeriesYlim,'k','visible','off');

% print out slice labels
xmax = min(get(gData.hTseriesAxis,'XTick'));
ymax = max(get(gData.hTseriesAxis,'yTick'));
for slicenum = 0:gData.spikeInfo.dim(3)
  if slicenum==0
    htext = text(xmax,ymax,'Slices:','parent',gData.hTseriesAxis);
  else
    htext = text(xmax,ymax,sprintf('%i',slicenum),'parent',gData.hTseriesAxis);
    set(htext,'Color',gData.colorOrder(slicenum,:));
  end
  position = get(htext,'Position');
  textextent = get(htext,'Extent');
  position(2) = position(2)-textextent(4);
  position(1) = position(1)+10;
  set(htext,'Position',position);
  xmax = xmax+textextent(3);
end
ylabel(gData.hTseriesAxis,{'Mean over each slice','(% signal change)'});

% plot Z value Matrix and Spikes
gData.hSpikeMatrix = subplot(gData.nrows,gData.ncols,plot2,'parent',hFigure);
title(gData.hSpikeMatrix,'Spike matrix: Use mouse or arrow keys to move the cursor and display corresponding image, FFT and spikes');
colormap(gData.hSpikeMatrix,'gray');

guidata(hFigure,gData);
setNewScan(hFigure)

% display all spikes
if fieldIsNotDefined(gData.spikeInfo,'currentCriterion')
  gData.spikeInfo.currentCriterion=gData.spikeInfo.criterion;
end
setCriterion(hFigure,gData.spikeInfo.currentCriterion);

set(hFigure,'WindowButtonDownFcn',{@MouseDownCallback});
set(hFigure,'WindowButtonMotionFcn',{@MouseMotionCallback});
set(hFigure,'KeyPressFcn',{@KeypressCallback});
set(hFigure,'interruptible','off')
set(hFigure,'BusyAction','cancel');

%%%%%%%%%%%%%%%%%%%%
%%   setNewScan   %%
%%%%%%%%%%%%%%%%%%%%
function setNewScan(hFigure)

gData = guidata(hFigure);

%change time-series Y scale
set(gData.hTseriesAxis,'xLim',[0 gData.spikeInfo.dim(4)]);
gData.tSeriesYlim = [floor(min(gData.spikeInfo.sliceMeans(:))) ceil(max(gData.spikeInfo.sliceMeans(:)))];
set(gData.hTseriesAxis,'yLim',gData.tSeriesYlim);
set(gData.hCursorTseries,'yData',gData.tSeriesYlim,'visible','on');

% plot Z value Matrix
hold(gData.hSpikeMatrix,'off');
if isfield(gData.spikeInfo,'maxZvalue')
  imagesc(gData.spikeInfo.maxZvalue,'parent',gData.hSpikeMatrix);
else
  imagesc(zeros(gData.spikeInfo.dim(3),gData.spikeInfo.dim(4)),'parent',gData.hSpikeMatrix);
end
hold(gData.hSpikeMatrix,'on');
caxis(gData.hSpikeMatrix,[0 8]);
set(gData.hSpikeMatrix,'xLim',[0.5 gData.spikeInfo.dim(4)+.5],'yLim',[0.5 gData.spikeInfo.dim(3)+.5]);
set(gData.hSpikeMatrix,'ydir','normal');
xlabel(gData.hSpikeMatrix,sprintf('Frame number'));
ylabel(gData.hSpikeMatrix,'Slice number');
h = colorbar('peer',gData.hSpikeMatrix,'location','southoutside');
set(get(h,'XLabel'),'string','Max FFT Z value in Slice/Frame')


if ~gData.spikeInfo.n
    gData.spikeInfo.currentSpike =0;
elseif fieldIsNotDefined(gData.spikeInfo,'currentSpike')
    gData.spikeInfo.currentSpike =1;
end
if fieldIsNotDefined(gData.spikeInfo,'currentFrame') || fieldIsNotDefined(gData.spikeInfo,'currentSlice')
  if gData.spikeInfo.currentSpike
    gData.spikeInfo.currentFrame =gData.spikeInfo.time(gData.spikeInfo.currentSpike);
    gData.spikeInfo.currentSlice =gData.spikeInfo.slice(gData.spikeInfo.currentSpike);
  else
    gData.spikeInfo.currentFrame =round(gData.spikeInfo.dim(4)/2);
    gData.spikeInfo.currentSlice =round(gData.spikeInfo.dim(3)/2);
  end    
end
if fieldIsNotDefined(gData.spikeInfo,'currentCriterion')
  gData.spikeInfo.currentCriterion=gData.spikeInfo.criterion;
else
  gData.spikeInfo.currentCriterion=max(gData.spikeInfo.criterion,gData.spikeInfo.currentCriterion);
end

guidata(hFigure,gData);
setCriterion(hFigure,gData.spikeInfo.currentCriterion);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   MouseDownCallback     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MouseDownCallback(hFigure,eventData)

gData = guidata(hFigure);

spikeMatrixCoords = round(get(gData.hSpikeMatrix,'CurrentPoint'));
if all(spikeMatrixCoords(1,[1 2])>0 & spikeMatrixCoords(1,[1 2])<=gData.spikeInfo.dim([4 3]))
  gData.spikeInfo.currentSlice = spikeMatrixCoords(1,2);
  gData.spikeInfo.currentFrame = spikeMatrixCoords(1,1);
  guidata(hFigure,gData);
  spikePlotImage(hFigure,0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   MouseMotionCallback     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MouseMotionCallback(hFigure,eventData)

gData = guidata(hFigure);

spikeMatrixCoords = round(get(gData.hSpikeMatrix,'CurrentPoint'));
tSeriesCoords = round(get(gData.hTseriesAxis,'CurrentPoint'));
% if ishandle(gData.hCursorTseries)
%   delete(gData.hCursorTseries);
% end
if all(spikeMatrixCoords(1,[1 2])>0 & spikeMatrixCoords(1,[1 2])<=gData.spikeInfo.dim([4 3])) ||...
  all(tSeriesCoords(1,[1 2])>[0 gData.tSeriesYlim(1)] & tSeriesCoords(1,[1 2])<=[gData.spikeInfo.dim(4) gData.tSeriesYlim(2)])
  set(gData.hCursorTseries,'xData',repmat(tSeriesCoords(1,1),1,2));
end

guidata(hFigure,gData);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   KeypressCallback     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function KeypressCallback(hFigure,eventData)

gData = guidata(hFigure);
pressedKey=double(get(hFigure,'CurrentCharacter'));
if ~isempty(pressedKey)
  switch(pressedKey)
    case 28 %left arrow
      gData.spikeInfo.currentFrame = max(gData.spikeInfo.currentFrame-1,1);
    case 29 %right arrow
      gData.spikeInfo.currentFrame = min(gData.spikeInfo.currentFrame+1,gData.spikeInfo.dim(4));
    case 30 %up arrow
      gData.spikeInfo.currentSlice = min(gData.spikeInfo.currentSlice+1,gData.spikeInfo.dim(3));
    case 31 %down arrow
      gData.spikeInfo.currentSlice = max(gData.spikeInfo.currentSlice-1,1);
    otherwise
      return;
  end
  guidata(hFigure,gData);
  spikePlotImage(hFigure,0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   spikePlotController   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spikePlotController(hFigure)

gData = guidata(hFigure);

% check all scans that have spikeInfo
spikeInfoScans = {};
for i = 1:viewGet(gData.v,'nScans',gData.spikeInfo.groupNum);
  if i==gData.spikeInfo.scanNum %for the current scan, get the current spikeInfo, not the one from the view
    thisSpikeinfo=gData.spikeInfo;
  else
    thisSpikeinfo = viewGet(gData.v,'spikeInfo',i,gData.spikeInfo.groupNum);
  end
  if ~isempty(thisSpikeinfo);
    if fieldIsNotDefined(thisSpikeinfo,'currentCriterion')
      criterion=thisSpikeinfo.criterion;
      spikeNumber = thisSpikeinfo.n;
    else
      criterion=thisSpikeinfo.currentCriterion;
      spikeNumber = thisSpikeinfo.currentSpikeNumber;
    end
    spikeInfoScans{end+1} = sprintf('%i: %s (%i spikes at z>%.2f)',...
      i,viewGet(gData.v,'description',i,thisSpikeinfo.groupNum),...
      spikeNumber,criterion);
    if i==gData.spikeInfo.scanNum
      thisSpikeInfoScans = spikeInfoScans{end};
    end
  end
end
spikeInfoScans = putOnTopOfList(thisSpikeInfoScans,spikeInfoScans);

% now put up a control dialog
paramsInfo{1}  = {'scanNum',spikeInfoScans,'Scan number to view'};
paramsInfo{end+1} = {'recompute',[],'type=pushbutton','buttonString=Recompute spike detection','callback',@spikePlotRecomputeCallback,'callbackArg',hFigure,'Recomputer spike detection'};
% if gData.spikeInfo.n > 0
%   paramsInfo{end+1} = {'spikeNum',1,sprintf('minmax=[1 %i]',gData.spikeInfo.n),'incdec=[-1 1]','round=1','Which spike to display'};
paramsInfo{end+1} = {'criterion',max(gData.spikeInfo.criterion,gData.spikeInfo.currentCriterion),sprintf('minmax=[%i inf]',gData.spikeInfo.criterion),'incdec=[-.5 .5]','To change the value of the criterion used (must be > criterion used for computing the detection)'};
% else
%   paramsInfo{end+1} = {'noSpikes','No spikes found','editable=0','type=string','No spikes found in scan'};
% end  
mrParamsDialog(paramsInfo,'mlrSpikeDetector',[],@spikePlotCallback,hFigure,{@spikePlotOKCallback,hFigure});


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   spikePlotCallback   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spikePlotCallback(params,hFigure)

gData = guidata(hFigure);

% get the scanNum
scanNum = str2num(strtok(params.scanNum,':'));

if scanNum ~= gData.spikeInfo.scanNum
  %save the current gData.spikeInfo
  gData.v = viewSet(gData.v,'spikeInfo',gData.spikeInfo,gData.spikeInfo.scanNum,gData.spikeInfo.groupNum);
  % load up the new spike info
  gData.spikeInfo = viewGet(gData.v,'spikeInfo',scanNum,gData.spikeInfo.groupNum);
  guidata(hFigure,gData);
  setNewScan(hFigure);
  spikePlotController(hFigure);
else
  if isfield(params,'criterion') && params.criterion~=gData.spikeInfo.currentCriterion
    setCriterion(hFigure,params.criterion);
    %change scan description
    %relaunch spikePlotcontroller (cannot use mrParamsSet to change the string of a popupmenu...)
    spikePlotController(hFigure);
%   elseif isfield(params,'spikeNum') && params.spikeNum~=gData.spikeInfo.currentSpike
%     spikePlotImage(hFigure,params.spikeNum);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   spikePlotOKCallback   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spikePlotOKCallback(hFigure)

gData = guidata(hFigure);

%save the current spikeInfo
gData.v = viewSet(gData.v,'spikeInfo',gData.spikeInfo,gData.spikeInfo.scanNum,gData.spikeInfo.groupNum);
%here there might be a problem if the view has changed outside mlrSpikeDetector...
closeGraphWin(hFigure);


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   spikePlotCallback   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = spikePlotRecomputeCallback(hFigure)

val = [];
gData = guidata(hFigure);
%save the current spikeInfo
gData.v = viewSet(gData.v,'spikeInfo',gData.spikeInfo,gData.spikeInfo.scanNum,gData.spikeInfo.groupNum);
closeGraphWin(hFigure);

mlrSpikeDetector(gData.v,gData.spikeInfo.scanNum,gData.spikeInfo.groupNum,'recompute=1');


%%%%%%%%%%%%%%%%%%%%%%%%
%%   setCriterion   %%
%%%%%%%%%%%%%%%%%%%%%%%%
function setCriterion(hFigure,thisCriterion)

gData = guidata(hFigure);

gData.spikeInfo.currentCriterion = thisCriterion;

if isfield(gData,'hBox') && all(ishandle(gData.hBox(:)))
  delete(gData.hBox);
  gData = rmfield(gData,'hBox');
end

%draw a contour around each slice/frame that contains spike over the current criterion
boxXcoords = [-.5 .5;-.5 .5;-.5 -.5;.5 .5];%;-.5 .5];
boxYcoords = [-.5 -.5;.5 .5;-.5 .5;-.5 .5];%;.5 -.5];
cSpike=0;
if gData.spikeInfo.n && isfield(gData.spikeInfo.spikelocs{1},'zValue')
  for iSpike = 1:gData.spikeInfo.n
    if any(gData.spikeInfo.spikelocs{iSpike}.zValue>thisCriterion)
      cSpike = cSpike+1;
      gData.hBox(cSpike,:) = plot(gData.hSpikeMatrix,gData.spikeInfo.time(iSpike)+boxXcoords,gData.spikeInfo.slice(iSpike)+boxYcoords,'g');
    end
  end
else  %for an old spikeInfo structure, we don't have the zValue info
  mrWarnDlg('(mlrSpikeDetector:setCriterion) There is no zValue information in spikeInfo, please recompute spike detection for this scan.');
  for cSpike = 1:gData.spikeInfo.n
      gData.hBox(cSpike,:) = plot(gData.hSpikeMatrix,gData.spikeInfo.time(cSpike)+boxXcoords,gData.spikeInfo.slice(cSpike)+boxYcoords,'g');
  end
end

gData.spikeInfo.currentSpikeNumber = cSpike;
title(gData.hTseriesAxis,...
      {sprintf('%i spikes found at z>%0.1f)',cSpike,thisCriterion),...
      '(Note that spike detection is done on each FFT component not on the mean of the slice)'},...
      'interpreter','none');

guidata(hFigure,gData); 

% display images with possible artifacts for the current slice/frame
spikePlotImage(hFigure,gData.spikeInfo.currentSpike);

%%%%%%%%%%%%%%%%%%%%%%%%
%%   spikePlotImage   %%
%%%%%%%%%%%%%%%%%%%%%%%%
function spikePlotImage(hFigure,spikeNum)

gData = guidata(hFigure);

%if there is no spikenum passed (0), see if there are spikes for this slice and frame
if ~spikeNum
  spikeNum = find(gData.spikeInfo.slice==gData.spikeInfo.currentSlice & gData.spikeInfo.time==gData.spikeInfo.currentFrame);
  if isempty(spikeNum)
    spikeNum=0;
  end
else
%otherwise, set the the current slice and frame on this spike
  gData.spikeInfo.currentFrame =gData.spikeInfo.time(spikeNum);
  gData.spikeInfo.currentSlice =gData.spikeInfo.slice(spikeNum);
end
gData.spikeInfo.currentSpike = spikeNum;
thisFrame = gData.spikeInfo.currentFrame;
thisSlice = gData.spikeInfo.currentSlice;

%draw the cursor
if isfield(gData,'hCursor') && all(ishandle(gData.hCursor))
  delete(gData.hCursor);
end
boxXcoords = [-.5 .5;-.5 .5;-.5 -.5;.5 .5];%;-.5 .5];
boxYcoords = [-.5 -.5;.5 .5;-.5 .5;-.5 .5];%;.5 -.5];
gData.hCursor = plot(gData.hSpikeMatrix,thisFrame+boxXcoords,thisSlice+boxYcoords,'m','linewidth',3);
%guidata(hFigure,gData); %save the new handle immediately
drawnow

%draw the mean time-series
if isfield(gData,'hTseries') && all(ishandle(gData.hTseries))
  delete(gData.hTseries);
end
for slicenum = [thisSlice+1:gData.spikeInfo.dim(3) 1:thisSlice-1]
  gData.hTseries(slicenum) = plot(gData.hTseriesAxis,gData.spikeInfo.sliceMeans(slicenum,:),'.-','color',gData.colorOrder(slicenum,:));
end
gData.hTseries(thisSlice) = plot(gData.hTseriesAxis,gData.spikeInfo.sliceMeans(thisSlice,:),'.-','color',gData.colorOrder(thisSlice,:),'lineWidth',2);

% display image
% read data
gData.v = viewSet(gData.v,'curGroup',gData.spikeInfo.groupNum);
data = loadTSeries(gData.v,gData.spikeInfo.scanNum,thisSlice,thisFrame);
imageg(data,0.5,gData.hImagePlot);
% and title
title(gData.hImagePlot,{sprintf('Slice %i, Frame %i',thisSlice,thisFrame),'(Normalized intensity values)'});
data(isnan(data)) = nanmean(data(:));
fftimage = abs(fftshift(fft2(data)));
% set gamma
imageg(fftimage,0.5,gData.hFftPlot);
title(gData.hFftPlot,{'FFT',sprintf('Components highlighted in green: z>%.2f',gData.spikeInfo.currentCriterion)});

if spikeNum
  hold(gData.hFftPlot,'on')
  % draw boxes around noise values
  if isfield(gData.spikeInfo.spikelocs{1},'zValue')
    for iLoc = 1:length(gData.spikeInfo.spikelocs{spikeNum}.x)
      if gData.spikeInfo.spikelocs{spikeNum}.zValue(iLoc)>gData.spikeInfo.currentCriterion
        plot(gData.hFftPlot,gData.spikeInfo.spikelocs{spikeNum}.y(iLoc)+boxXcoords,gData.spikeInfo.spikelocs{spikeNum}.x(iLoc)+boxYcoords,'g','linewidth',1); %why should x and y be inverted (and does that matter ?)
      end
    end
  else  %for an old spikeInfo structure, we don't have the zValue info
    mrWarnDlg('(mlrSpikeDetector:setCriterion) There is no zValue information in spikeInfo, please recompute spike detection for this scan.');
    for iLoc = 1:length(gData.spikeInfo.spikelocs{spikeNum}.x)
      plot(gData.hFftPlot,gData.spikeInfo.spikelocs{spikeNum}.y(iLoc)+boxXcoords,gData.spikeInfo.spikelocs{spikeNum}.x(iLoc)+boxYcoords,'g','linewidth',1); %why should x and y be inverted (and does that matter ?)
    end
  end
  hold(gData.hFftPlot,'off')
end

guidata(hFigure,gData);


% imageg.m
%
%      usage: imageg(data array,gamma)
%      usage: imageg(data structure,slicenum,volumenum,gamma)
%         by: justin gardner
%       date: 05/09/03
%    purpose: gamma correct image display for epi images
%
function d = imageg(d,gamma,handle)

if size(d,2) < size(d,1)
  d = d';
end

% scale image values to between 0 and 1
imagemax = max(max(d));
imagemin = min(min(d));
d = (d-imagemin)./(imagemax-imagemin);

% apply monitor gamma
d = d.^gamma;

% rescale to 0-255 uint8
d = floor(255*d);

image(d,'parent',handle);
colormap(handle,gray(256));
axis(handle,'off');
axis(handle,'equal');
