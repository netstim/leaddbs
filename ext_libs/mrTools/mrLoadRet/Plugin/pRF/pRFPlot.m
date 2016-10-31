% pRFPlot.m
%
%        $Id:$ 
%      usage: pRFPlot(v,overlayNum,scanNum,x,y,z,,roi)
%         by: justin gardner
%       date: 11/22/11
%    purpose: plot function for displaying results of pRF analysis
%
function pRFPlot(v,overlayNum,scanNum,x,y,z,roi)

% check arguments
if ~any(nargin == [7])
  help pRFPlot
  return
end

% see if the shift key is down
%shiftDown = any(strcmp(get(viewGet(v,'figureNumber'),'CurrentModifier'),'shift'));
shiftDown = any(strcmp(get(viewGet(v,'figureNumber'),'SelectionType'),'extend'));

% check if pRF has been run
a = viewGet(v,'Analysis');
if ~isfield(a,'type') || ~strcmp(a.type,'pRFAnal')
  disp(sprintf('(pRFPlot) pRF analysis has not been run on this scan'));
  return
end

% get the d
d = viewGet(v,'d',scanNum);
if isempty(d),disp(sprintf('(pRFPlot) Could not find d structure for this scan'));return,end

% get the parametrs of the pRF fit
r2 = viewGet(v,'overlayData',scanNum,viewGet(v,'overlayNum','r2'));
if isempty(r2)
  disp(sprintf('(pRFPlot) pRF analysis has not been run on this scan'));
  return
end
thisR2 = r2(x,y,z);
polarAngle = viewGet(v,'overlayData',scanNum,viewGet(v,'overlayNum','polarAngle'));
thisPolarAngle = polarAngle(x,y,z);
eccentricity = viewGet(v,'overlayData',scanNum,viewGet(v,'overlayNum','eccentricity'));
thisEccentricity = eccentricity(x,y,z);
rfHalfWidth = viewGet(v,'overlayData',scanNum,viewGet(v,'overlayNum','rfHalfWidth'));
thisRfHalfWidth = rfHalfWidth(x,y,z);

% roi
if ~shiftDown
  pRFPlotROI(v,roi,d,a,r2,eccentricity,polarAngle,rfHalfWidth);
end

% get the params that have been run
scanDims = viewGet(v,'scanDims',scanNum);
whichVoxel = find(d.linearCoords == sub2ind(scanDims,x,y,z));
r = d.r(whichVoxel,:);

% if no voxel has been found in precomputed analysis then do fit (or if shift is down)
if isempty(whichVoxel) || shiftDown
  % check if shift is being held down, in which case we reget parameters
  if shiftDown
    fit = pRFFit(v,overlayNum,scanNum,x,y,z,roi);
  else
    fit = pRFFit(v,overlayNum,scanNum,x,y,z,roi,'fitTypeParams',a.params.pRFFit);
  end
  if isempty(fit),return,end
  % set the overlays
  r2(x,y,z) = fit.r2;
  polarAngle(x,y,z) = fit.polarAngle;
  eccentricity(x,y,z) = fit.eccentricity;
  rfHalfWidth(x,y,z) = fit.std;
  % reset the overlays
  v = viewSet(v,'overlayDataReplace',r2,'r2');
  v = viewSet(v,'overlayDataReplace',polarAngle,'polarAngle');
  v = viewSet(v,'overlayDataReplace',eccentricity,'eccentricity');
  v = viewSet(v,'overlayDataReplace',rfHalfWidth,'rfHalfWidth');
  % now refresh the display
  refreshMLRDisplay(viewGet(v,'viewNum'));
  return
end

params = d.params(:,whichVoxel);
if isfield(d,'paramsInfo')
  paramsInfo = d.paramsInfo;
else
  paramsInfo = [];
end

% get params
m = pRFFit(v,scanNum,x,y,z,'stim',d.stim,'getModelResponse=1','params',params,'concatInfo',d.concatInfo,'fitTypeParams',a.params.pRFFit,'paramsInfo',paramsInfo);
% and plot, set a global so that we can use the mouse to display
% different time points
global gpRFPlot;
gpRFPlot.fignum = selectGraphWin;

% clear callbacks
set(gpRFPlot.fignum,'WindowButtonMotionFcn','');

% keep the stim
gpRFPlot.d = d;
gpRFPlot.rfModel = m.rfModel;

% keep the axis that has the time series
gpRFPlot.a = subplot(5,5,[1:4 6:9 11:14 16:19]);
% plot the rectangle that shows the current stimuli
% FIX: Start time
gpRFPlot.t = 50;
gpRFPlot.hRect = rectangle('Position',[gpRFPlot.t-4 min(m.tSeries) 4 max(m.tSeries)-min(m.tSeries)],'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]);
hold on
% plot time series
plot(m.tSeries,'k.-');
axis tight
% plot model
plot(m.modelResponse,'r-');
if d.concatInfo.n > 1
  vline(d.concatInfo.runTransition(2:end,1));
end
xlabel('Time (volumes)');
ylabel('BOLD (%)');
% convert coordinates back to x,y for display
[thisx thisy] = pol2cart(thisPolarAngle,thisEccentricity);
title(sprintf('[%i %i %i] r^2=%0.2f polarAngle=%0.2f eccentricity=%0.2f rfHalfWidth=%0.2f %s [x=%0.2f y=%0.2f]\n%s',x,y,z,thisR2,r2d(thisPolarAngle),thisEccentricity,thisRfHalfWidth,a.params.pRFFit.rfType,thisx,thisy,num2str(r,'%0.2f ')));
% plot the rf
a = subplot(5,5,[10 15 20]);
imagesc(d.stimX(:,1),d.stimY(1,:),flipud(m.rfModel'));
set(a,'Box','off');
set(a,'Color',[0.8 0.8 0.8]);
set(a,'TickDir','out');
axis equal
axis tight
hold on
hline(0,'w:');vline(0,'w:');
% plot the canonical
subplot(5,5,5);cla
plot(m.canonical.time,m.canonical.hrf,'k-');
title(sprintf('lag: %0.2f tau: %0.2f',m.p.canonical.timelag,m.p.canonical.tau));

% display the stimulus images
plotStim(gpRFPlot.t);

% now set callback
set(gpRFPlot.fignum,'WindowButtonMotionFcn',@pRFPlotMoveMouse);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    mlrprFPloMoveMouse    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pRFPlotMoveMouse(hWindow,event)

global gpRFPlot;
if ~ishandle(gpRFPlot.a),return,end

currentPoint = get(gpRFPlot.a ,'CurrentPoint');
coord = round(currentPoint(1,1:2));
a = axis(gpRFPlot.a);
if (coord(1) >= a(1)) && (coord(1) <= a(2)) && (coord(2) >= a(3)) && (coord(2) <= a(4))
  % move rectangle
  pos = get(gpRFPlot.hRect,'Position');
  pos(1) = coord(1)-4;
  set(gpRFPlot.hRect,'Position',pos);
  % redisplay stimulus images
  plotStim(coord(1))
end

%%%%%%%%%%%%%%%%%%
%    plotStim    %
%%%%%%%%%%%%%%%%%%
function plotStim(t)

global gpRFPlot;

for i = 1:5
  a = subplot(5,5,20+i,'Parent',gpRFPlot.fignum);
  cla(a);
  thist = t-5+i;
  if thist >= 1 
    im = [];
    % get the scan and volume
    thisScan = gpRFPlot.d.concatInfo.whichScan(thist);
    thisVolume = gpRFPlot.d.concatInfo.whichVolume(thist);
    junkFrames = gpRFPlot.d.concatInfo.totalJunkedFrames(thisScan);
    im(:,:,3) = flipud(0.7*gpRFPlot.d.stim{thisScan}.im(:,:,thisVolume+junkFrames)');
    im(:,:,2) = flipud(0.7*gpRFPlot.d.stim{thisScan}.im(:,:,thisVolume+junkFrames)');
    im(:,:,1) = flipud(0.7*gpRFPlot.d.stim{thisScan}.im(:,:,thisVolume+junkFrames)'+0.3*gpRFPlot.rfModel');
    % swap and flip so that it will display correctly
    image(gpRFPlot.d.stimX(:,1),gpRFPlot.d.stimY(1,:),im,'Parent',a);
    axis image
    hold(a,'on');
    hline(0,'w:',a);vline(0,'w:',a);
    title(a,sprintf('t=%i',thist));
  end
end

%%%%%%%%%%%%%%%%%%%%
%    pRFPlotROI    %
%%%%%%%%%%%%%%%%%%%%
function pRFPlotROI(v,roi,d,a,r2,eccentricity,polarAngle,rfHalfWidth)

if length(roi)
  % check for already plotted
  minr2 = viewGet(v,'overlayMin','r2');
  scanNum = viewGet(v,'curScan');
  groupNum = viewGet(v,'curGroup');
  global gpRFPlotROI
  checkParams = {'roi','minr2','a','scanNum','groupNum'};
  replot = false;
  % if shift key is down then replot
  f = viewGet(v,'fignum');
  if ~isempty(f) && any(strcmp(get(f,'CurrentModifier'),'shift')),replot=true;end
  for i = 1:length(checkParams)
    if ~isfield(gpRFPlotROI,checkParams{i}) || ~isequal(gpRFPlotROI.(checkParams{i}),eval(checkParams{i}))
      replot = true;
    end
    gpRFPlotROI.(checkParams{i}) = eval(checkParams{i});
  end
  if ~replot, return, end
  disp(sprintf('(pRFPlot) Displaying ROI fig'));
  mlrSmartfig('pRFPlotROI','reuse');clf
  
  minX = min(d.stimX(:));
  maxX = max(d.stimX(:));
  minY = min(d.stimY(:));
  maxY = max(d.stimY(:));
  
  % see what kind of fit we have.
  if strcmp(a.params.pRFFit.rfType,'gaussian-hdr')
    % plot also the hdr parameters
    numRowsPerROI = 2;
    numCols = 3;
    % set up fields for plotting extra hdr parameters
    if a.params.pRFFit.diffOfGamma
      plotParams = [4 5 6 7 8];
      plotParamsNames = {'timelag','tau','amplitudeRatio','timelag2','tau2'};
      numCols = 5;
    else
      plotParams = [4 5];
      plotParamsNames = {'timelag','tau'};
    end
  else
    numRowsPerROI = 1;
    numCols = 3;
    plotParams = [];
    plotParamsNames = {};
  end
  
  
  for roiNum = 1:length(roi)
    % get coordinates
    % roiCoords = getROICoordinates(v,roi{roiNum},[],[],'straightXform=1');
    roiCoords = getROICoordinates(v,roi{roiNum});
    roiCoordsLinear = sub2ind(viewGet(v,'scanDims'),roiCoords(1,:),roiCoords(2,:),roiCoords(3,:));
    % get values for the roi
    thisr2 = r2(roiCoordsLinear);
    % only use voxels above current r2 min
    roiCoordsLinear = roiCoordsLinear(find(thisr2 >minr2));
    % sort them
    [thisr2sorted r2index] = sort(r2(roiCoordsLinear));
    roiCoordsLinear = roiCoordsLinear(r2index);
    % get values for these voxels
    thisr2 = r2(roiCoordsLinear);
    thisEccentricity = eccentricity(roiCoordsLinear);
    thisPolarAngle = polarAngle(roiCoordsLinear);
    thisRfHalfWidth = rfHalfWidth(roiCoordsLinear);
    % convert to cartesian
    [thisX thisY] = pol2cart(thisPolarAngle,thisEccentricity);
    c = [1 1 1];
    
    % plot RF coverage
    subplot(length(roi)*numRowsPerROI,numCols,1+(roiNum-1)*numCols*numRowsPerROI);
    for i = 1:length(thisX)
      if ~isnan(thisr2(i))
        plotCircle(thisX(i),thisY(i),thisRfHalfWidth(i),1-c*thisr2(i)/max(thisr2));
        hold on
      end
    end
    xaxis(minX,maxX);
    yaxis(minY,maxY);
    axis square
    hline(0);
    vline(0);
    xlabel('x (deg)');
    ylabel('y (deg)');
    title(sprintf('%s rf (r2 cutoff: %0.2f)',roi{roiNum}.name,minr2));
    
    % plot RF centers
    subplot(length(roi)*numRowsPerROI,numCols,2+(roiNum-1)*numCols*numRowsPerROI);
    for i = 1:length(thisX)
      if ~isnan(thisr2(i))
        plot(thisX(i),thisY(i),'k.','Color',1-c*thisr2(i)/max(thisr2), 'markersize', 10);
        hold on
      end
    end
    xaxis(minX,maxX);
    yaxis(minY,maxY);
    axis square
    hline(0);
    vline(0);
    xlabel('x (deg)');
    ylabel('y (deg)');
    title(sprintf('%s centers',roi{roiNum}.name));
    
    % plot eccentricity vs. rfHalfWidth
    subplot(length(roi)*numRowsPerROI,numCols,3+(roiNum-1)*numCols*numRowsPerROI);
    for i = 1:length(thisX)
      if ~isnan(thisr2(i))
        plot(thisEccentricity(i),thisRfHalfWidth(i),'k.','Color',1-c*thisr2(i)/max(thisr2), 'markersize', 10);
        hold on
      end
    end
    hold on
    % limit the fit to the central 6 deg (b/c it is often off for higher eccentricities)
    eccLimit = 6;
    ind = thisEccentricity <= eccLimit;
    if any(ind)
%      regfit = myregress(thisEccentricity(ind),thisRfHalfWidth(ind),0,0);
      w = diag(thisr2(ind));
      x = thisEccentricity(ind);
      x = [x(:) ones(size(x(:)))];
      y = thisRfHalfWidth(ind);
      beta = ((x'*w*x)^-1)*(x'*w)*y';
      maxXaxis = min(maxX,maxY);
      xaxis(0,maxXaxis);
      yaxis(0,maxXaxis);
      if ~isempty(beta)
	plot([0 maxXaxis],[0 maxXaxis]*beta(1)+beta(2),'k-');
      end
      xlabel('Eccentricity (deg)');
      ylabel('RF half width (deg)');
%      title(sprintf('slope: %0.2f (%s) offset: %0.2f (%s) (r2=%0.2f)',beta(1),pvaldisp(regfit.pm),beta(2),pvaldisp(regfit.pb),regfit.r2));
      axis square
    else
      disp(sprintf('(pRFPlot) No matching fits to plot with eccentricity less than %f',eccLimit));
    end
    % plot hdr parameters, first get the voxels to plot
    [temp dCoords] = intersect(d.linearCoords,roiCoordsLinear);
    for i = 1:length(plotParams)
      subplot(length(roi)*numRowsPerROI,numCols,numCols+i+(roiNum-1)*numCols*numRowsPerROI);
      hist(d.params(plotParams(i),dCoords));
      xlabel(plotParamsNames{i});
      ylabel('n');
      if exist('plotmean')==2
        plotmean(d.params(plotParams(i),dCoords));
      end
    end
  end
end

%%%%%%%%%%%%%%%%%%%%
%    plotCircle    %
%%%%%%%%%%%%%%%%%%%%
function h = plotCircle(xCenter,yCenter,radius,c)

a = 0:0.01:2*pi;
h = plot(xCenter+radius*cos(a),yCenter+radius*sin(a),'k-','Color',c);


%%%%%%%%%%%%%
%%   r2d   %%
%%%%%%%%%%%%%
function degrees = r2d(angle)

degrees = (angle/(2*pi))*360;

% if larger than 360 degrees then subtract
% 360 degrees
while (sum(degrees>360))
  degrees = degrees - (degrees>360)*360;
end

% if less than 360 degreees then add 
% 360 degrees
while (sum(degrees<-360))
  degrees = degrees + (degrees<-360)*360;
end

