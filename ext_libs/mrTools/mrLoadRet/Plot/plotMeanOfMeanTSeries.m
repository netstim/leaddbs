function plotMeanOfMeanTSeries(view, groupNum, roiList, scanList, varargin)
%
% plotMeanOfMeanTSeries(view, groupNum, roiList, scanList, [param1], [value1], [param2], [value2])
%
% Computes mean tseries for each ROI and each scan, and then plots the mean
% and SEM across scans (i.e., one subplot for each ROI).
%
% groupNum: group number (default: current group)
% roiList: vector of ROI numbers (default: [1:nROIs])
% scanList: vector of scan numbers (default: [1:nscans])
% other arguments as in percentTSeries
%
% djh 9/2005

if ieNotDefined('groupNum')
	groupNum = viewGet(view,'currentGroup');
end
view = viewSet(view,'currentGroup',groupNum);
if ieNotDefined('roiNUM')
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

% Compute mean tseries
tseries = meanTSeries(view, groupNum, roiList, scanList, varargin{:});

% Select window
selectGraphWin;
fontSize = 14;
set(gca,'FontSize',fontSize);

tseriesROI = zeros(nscans,nframes);
for iROI = 1:nROIs
	roiName = viewGet(view,'roiName',roiList(iROI));
	for iscan = 1:nscans
		scan = scanList(iscan);
		if (nframes ~= viewGet(view,'nFrames',scan,groupNum))
			mrErrorDlg('Cannot average these scans because they have different numFrames.');
		end
		if (framePeriod ~= viewGet(view,'framePeriod',scan,groupNum))
			mrWarnDlg('These scans  have different frame periods.');
		end
		tseriesROI(iscan,:) = tseries{iROI,iscan}';
	end
	tseriesMean = mean(tseriesROI);
	tseriesStd = std(tseriesROI,0);
	tseriesSEM = tseriesStd / sqrt(nscans);

	% Plot it
	subplot(nROIs,1,iROI);
	time = linspace(framePeriod,nframes*framePeriod,nframes);
	dsErrorsurface(time, tseriesMean, tseriesSEM, [0.5 0.5 0.5], 1);
	hold on
	plot(time,tseriesMean,'k','LineWidth',2);
	hold off

	% Ticks and labels
	headerStr = ['Times series from ',roiName,', group ',groupName,', scans: ',num2str(scanList)];
	title(headerStr,'interpreter','none');
	set(gca,'XLim',ceil([0,nframes*framePeriod]));
	xlabel('Time (sec)');
	ylabel('fMRI response');
end

return
