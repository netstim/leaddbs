function plotMeanTSeries(view, groupNum, roiList, scanList, varargin)
%
% plotMeanTSeries(view, groupNum, roiList, scanList, [param1], [value1], [param2], [value2])
% 
% Computes mean tseries for each ROI and each scan, and plots them.
%
% groupNum: default current group
% roiList: vector of ROI numbers (default: [1:nROIs])
% scanList: vector of scan numbers (default: [1:nscans])
% other arguments as in percentTSeries  
%
% djh 9/2005

if ieNotDefined('groupNum')
	groupNum = viewGet(view,'currentGroup');
end
view = viewSet(view,'currentGroup',groupNum);
if ieNotDefined('roiList')
	roiList = [1:viewGet(view,'numberofROIs')];
end
if ieNotDefined('scanList')
	scanList = [1:viewGet(view,'nscans')];
end
nROIs = length(roiList);
nscans = length(scanList);

% Compute mean tseries
tseries = meanTSeries(view, groupNum, roiList, scanList, varargin{:});

% Select window
selectGraphWin;
set(gcf,'name',['Mean timeseries, group ' viewGet(view,'groupName')]);

for iROI = 1:nROIs
  roiName = viewGet(view,'roiName',roiList(iROI));
  for iscan = 1:nscans
    scan = scanList(iscan);
    nframes = viewGet(view,'nFrames',scan);
    framePeriod = viewGet(view,'framePeriod',scan);
    ts = tseries{iROI,iscan};
		
    if ~isempty(ts)
      % Plot it
      subplot(nscans,nROIs,sub2ind([nROIs nscans],iROI,iscan));
      time = linspace(framePeriod,nframes*framePeriod,nframes)';
      plot(time,ts,'k.-');
      % Ticks and labels
      headerStr = [roiName,' Scan ',num2str(scan)];
      title(headerStr,'interpreter','none');
      set(gca,'XLim',ceil([0,nframes*framePeriod]));
      xlabel('Time (sec)');
      ylabel('fMRI response');
    end
  end
end

return
