function plotMeanFourierAmp(view, groupNum, roiList, scanList, varargin)
%
% plotMeanFourierAmp(view, groupNum, roiList, scanList, [param1], [value1], [param2], [value2])
%
% Computes Fourier Amplitude for each voxel in each ROI, for each scan, and
% then plots the mean across scans and voxels (i.e., one subplot for each
% ROI).
%
% groupNum: group number (default: current group)
% roiList: vector of ROI numbers (default: [1:nROIs])
% scanList: vector of scan numbers (default: [1:nscans])
% other arguments as in percentTSeries
%
% djh 9/2005
% jlg 11/2006 calculation of SNR

if ieNotDefined('groupNum')
	groupNum = viewGet(view,'currentGroup');
end
view = viewSet(view,'currentGroup',groupNum);
if ieNotDefined('roiList')
	roiList = viewGet(view,'currentROI');
end
if ieNotDefined('scanList')
	scanList = [1:viewGet(view,'nscans')];
end

nROIs = length(roiList);
nscans = length(scanList);
nframes = viewGet(view,'nFrames',scanList(1),groupNum);
framePeriod = viewGet(view,'framePeriod',scanList(1),groupNum);
groupName = viewGet(view,'groupName',groupNum);
currentAnalysis = viewGet(view,'currentAnalysis');
ncycles(1:max(scanList)) = 0;
if ~isempty(currentAnalysis)
  analysisParams = viewGet(view,'analysisParams',currentAnalysis);
  if isfield(analysisParams,'ncycles')
    ncycles(1:length(analysisParams.ncycles)) = analysisParams.ncycles;
  end
end

% Get the tseries for all voxels in each scan and each ROI
tseriesAll = tseriesROI(view, groupNum, roiList, scanList, varargin{:});

% this is an option so that we can average across all scans to
% make the meaurement, it needs an input argument, a gui item
% and should be tested before it can be used. -jlg
averageAcrossScans = 0;

% Loop through scans, computing fft of timeSeries and collecting
% them into one big matrix
maxy = 0;
selectGraphWin;
set(gcf,'name',['Mean Fourier amplitude, group ' viewGet(view,'groupName')]);

for iROI = 1:nROIs
  ffts = [];
  roiName = viewGet(view,'roiName',roiList(iROI));
  for iscan = 1:nscans
    scan = scanList(iscan);
    if averageAcrossScans
      if (nframes ~= viewGet(view,'nFrames',scan,groupNum))
	mrErrorDlg('Cannot average these scans because they have different numFrames.');
      end
      if (framePeriod ~= viewGet(view,'framePeriod',scan,groupNum))
	mrWarnDlg('These scans  have different frame periods.');
      end
      ffts = [ffts, abs(fft(tseriesAll{iROI,iscan}))];
    else
      ffts = abs(fft(tseriesAll{iROI,iscan}));
      selectGraphWin(1);
      subplot(nscans,nROIs,sub2ind([nROIs nscans],iROI,iscan));
      % plot and get maximum y scale so that we can rescale later
      maxy = max(maxy,dispFourierAmp(ffts,nframes,framePeriod,groupName,scanList(iscan),ncycles,roiName,viewGet(view,'description',scan,groupNum)));
    end
  end
  % if we are averaging acrross scans then display here
  if averageAcrossScans
    selectGraphWin(1);
    subplot(nROIs,1,iROI);
    maxy = max(maxy,dispFourierAmp(ffts,nframes,framePeriod,groupName,scanList,ncycles,roiName));
  end
end

% rescale all axis, so that they match
for iROI = 1:nROIs
  for iscan = 1:nscans
    selectGraphWin(1);
    if averageAcrossScans
      subplot(nROIs,1,iROI);
    else
      subplot(nscans,nROIs,sub2ind([nscans,nROIs],iscan,iROI));
    end
    % rescale axis
    yaxis(0,maxy);
  end
end

return

% Test
plotScansFourierAmp(MLR.views{1},[],[],[],'detrend','None','spatialNormalization','None')

function maxy = dispFourierAmp(ffts,nframes,framePeriod,groupName,scanNum,ncycles,roiName,description)

% Compute mean fourier amplitude
meanfft = nanmean(ffts,2);
meanfft = meanfft / (nframes/2);
frequencies = [0:nframes-1]/(nframes*framePeriod);
meanfft = meanfft(1:floor(nframes/2));
frequencies = frequencies(1:floor(nframes/2));

% get mean noise amplitude, by just taking average over top 1/3 of
% frequency values
noiseFrequencies = ceil(length(meanfft)*2/3):length(meanfft);
noiseAmplitude = meanfft(noiseFrequencies);
meanNoiseAmplitude = nanmean(noiseAmplitude);
lowestNoiseFrequency = frequencies(min(noiseFrequencies));
noiseFrequencies = frequencies(noiseFrequencies);

% get the signal amplitude
signalAmplitude = 0;signalFrequency = 0;
if ncycles(scanNum(1)) > 0
  signalFrequency = frequencies(ncycles(scanNum)+1);
  signalAmplitude = meanfft(ncycles(scanNum)+1);
end

% Plot it
fontSize = 14;
set(gca,'FontSize',fontSize);
plot(frequencies,meanfft,'k.-','LineWidth',1);
hold on
plot(noiseFrequencies,noiseAmplitude,'o','MarkerEdgeColor',[0.5 0.5 0.5]);
headerStr = sprintf('%s group: %s scan: %s roi: %s',description,groupName,num2str(scanNum),roiName);
% display SNR info only if signal frequecny is defined
if signalFrequency ~= 0
  plot(signalFrequency,signalAmplitude,'ro');
  headerStr = sprintf('%s\nSignal (%0.4f Hz) = %0.2f Noise (%0.4fHz and above) = %0.2f SNR=%0.3f',headerStr,signalFrequency,signalAmplitude,lowestNoiseFrequency,meanNoiseAmplitude,signalAmplitude/meanNoiseAmplitude);
else
  headerStr = sprintf('%s\nNoise (%0.4fHz and above) = %0.2f',headerStr,lowestNoiseFrequency,meanNoiseAmplitude);
end  
title(headerStr,'interpreter','none');
xlim([0 max(frequencies)]);
xlabel('Frequency (Hz)');
ylabel('Fourier Amplitude');
hold on
axis tight

a = axis;
maxy = a(4);
zoom on
drawnow