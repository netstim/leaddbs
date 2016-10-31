% timecoursePlot
%
%      usage: timecoursePlot()
%         by: justin gardner
%       date: 03/15/07
%    purpose: a default interrogator function that plots 
%             the time series and its fft
%
function timecoursePlot(view,overlayNum,scan,x,y,s,roi)


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

% check for nans in the timeseries
nanLocs = find(isnan(tSeries));
if ~isempty(nanLocs)
  disp(sprintf('(timecoursePlot) Found %i nans in locs: %s',length(nanLocs),mlrnum2str(nanLocs(:)','sigfigs=0','compact=1')));
  disp(sprintf('(timecoursePlot) Replacing nan points with mean of tSeries'));
  % set those locations to the mean so that we can show an approximate fft
  tSeriesNanMean = nanmean(tSeries);
  if ~isnan(tSeriesNanMean)
    tSeries(nanLocs) = tSeriesNanMean;
  else
    tSeries(nanLocs) = 0;
  end
end



% get the mean and trend
model = [(1:nFrames);ones(1,nFrames)]';
wgts = model \ tSeries;
fit = model*wgts;

subplot(2,1,1);
plot(tSeries, 'k.-');
hold on
plot(fit, '-', 'Color', [.5 .5 .5]);
title(sprintf('Voxel: [%i, %i, %i], mean=%0.2f, trend=%0.2f (%% sig change)',x,y,s,wgts(2), 100*wgts(1)/wgts(2)));
xlabel('Volumes');
ylabel('fMRI Signal');
axis tight;
% draw borders between rund
concatInfo = viewGet(view,'concatInfo',scan);
if ~isempty(concatInfo)
  vline(concatInfo.runTransition(2:end,1)-1,'r-');
end

subplot(2,1,2);
fftTSeries = fft(tSeries);
% set mean to zero
fftTSeries(1) = 0;
% plot it 
plot(1:(length(fftTSeries)/2)-1,abs(fftTSeries(2:length(fftTSeries)/2)),'k.-');
title(sprintf('Voxel: [%i %i %i]',x,y,s));
xlabel('FFT components');
ylabel('FFT of fMRI Signal');
axis tight;
zoom on

