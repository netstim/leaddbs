% timeSeriesStatsPlot
%
%      usage: timeSeriesStatsPlot()
%         by: justin gardner
%       date: 02/10/10
%    purpose: interrogator for timeSeriesStats
%
function timeSeriesStatsPlot(view,overlayNum,scan,x,y,s,roi)


% check arguments
if ~any(nargin == [1:7])
  help eventRelatedPlot
  return
end

% select the window to plot into
selectGraphWin;
global MLR;
fignum = MLR.graphFigure;

% turn off menu/title etc.
set(fignum,'NumberTitle','off');
set(fignum,'Name','timecoursePlot');

junkFrames = viewGet(view, 'junkFrames', scan);
nFrames = viewGet(view,'nFrames',scan);

tSeries = squeeze(loadTSeries(view,scan,s,[],x,y));
jSeries = tSeries(1:junkFrames);
tSeries = tSeries(junkFrames+1:junkFrames+nFrames);

nCols = 1;
if ~isempty(roi)
  % get the roi
  roi = roi{1};
  nCols = 2;
  % get time series for roi
  roi = loadROITSeries(view,roi)
  roi.tSeriesMean = mean(roi.tSeries);
  subplot(2,nCols,2);
  plot(roi.tSeriesMean,'k.-');
  title(sprintf('Mean tSeries of %s (n=%i)\nMean: %f Median: %f STD: %f mean/std: %f',roi.name,roi.n,mean(roi.tSeriesMean),median(roi.tSeriesMean),std(roi.tSeriesMean),mean(roi.tSeriesMean)/std(roi.tSeriesMean)));
  subplot(2,nCols,4);
  roi.fftTSeriesMean = fft(roi.tSeriesMean);
  roi.fftTSeriesMean(1) = 0;
  plot(abs(roi.fftTSeriesMean(1:length(roi.fftTSeriesMean)/2)),'k.-');
  title(sprintf('Mean FFT of %s',roi.name));
end

% get the overlays
titleStr = sprintf('Voxel: [%i, %i, %i]',x,y,s);
for iOverlay = 1:viewGet(view,'nOverlays');
  % get name
  overlayName = viewGet(view,'overlayName',iOverlay);
  % get value
  overlay = viewGet(view,'overlayData',scan,iOverlay);
  overlayValue = overlay(x,y,s);
  % add to titleStr
  titleStr = sprintf('%s %s=%f',titleStr,overlayName,overlayValue);
  if iOverlay == 3
    titleStr = sprintf('%s\n',titleStr);
  end
end

subplot(2,nCols,1);
plot(tSeries, 'k.-');
hold on

title(titleStr);
xlabel('Volumes');
ylabel('fMRI Signal');
axis tight;
% draw borders between rund
concatInfo = viewGet(view,'concatInfo',scan);
if ~isempty(concatInfo)
  vline(concatInfo.runTransition(2:end,1)-1,'r-');
end

subplot(2,nCols,nCols+1);
fftTSeries = fft(tSeries);
% set mean to zero
fftTSeries(1) = 0;
% plot it 
plot(abs(fftTSeries(1:length(fftTSeries)/2)),'k.-');
title(sprintf('Voxel: [%i %i %i]',x,y,s));
xlabel('FFT components');
ylabel('FFT of fMRI Signal');
axis tight;

