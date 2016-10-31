function view = timeSeriesStats(view,params)
%
% view = timeSeriesStats(view,[params])
% 
% Loops throughs scans, loads corresponding tSeries, and computes time
% series statistics for each voxel:
% - mean 
% - median 
% - standard deviation
% - max frame-to-frame difference
% - max difference from median
%
% params: Optional initial parameters. Default: user is prompted via
%    GUI. Params must be a a structure with all of the following fields. 
% params.groupName: group of scans that will be analyzed.
%    Default: current group of view.
% params.scanList: vector specifying which scans to compute.
%    Default: all of the scans.
%
%
% Examples:
%
% paramss.groupName = 'Raw';
% n = viewGet([],'nScans',1)
% params.scanList = [1:n];
% view = timeSeriesStats(view,params);
%
% view = timeSeriesStats(view);
%
%
% djh, 7/2007
% $Id$	

if ~isview(view)
    help timeSeriesStats
    mrErrorDlg('(timeSeriesStats) Invalid view.')
end

% Get analysis parameters from timeSeriesStatsGUI
if ieNotDefined('params')
  % Initialize analysis parameters with default values
  groupName = viewGet(view,'groupName');
  params.groupName = groupName;
  % find out all the scans in the current group that have more than 2 volumes
  % since we can't compute tSeriesStats for scans with less
  groupNum = viewGet(view,'groupNum',groupName);
  for iScan = 1:viewGet(view,'nScans',groupNum)
    nVols(iScan) = viewGet(view,'nVolumes',iScan);
    if nVols(iScan) <= 2
      disp(sprintf('(timeSeriesStats) Scan %s:%i has %i volumes - need at least 2 volumes to run tSeriesStats',groupName,iScan,nVols(iScan)));
    end
  end
  scanList = find(nVols>2);
  if length(scanList) > 0
    params.scanList = selectInList(view,'scans','Select scans for timeSeriesStats',scanList);
    if isempty(params.scanList),return,end
  else
    disp(sprintf('(timeSeriesStats) Could not find any scans to process'));
    return
  end
end

% Reconcile params with current status of group and ensure that params
% has the required fields.
params = defaultReconcileParams(params.groupName,params);

% Abort if params empty
if ieNotDefined('params')
  mrMsgBox('timeSeriesStats cancelled',1);
  return
end

% Change group
groupName = params.groupName;
curGroup = viewGet(view,'currentGroup');
groupNum = viewGet(view,'groupNum',groupName);
if (groupNum ~= curGroup)
	mrWarnDlg(['Changing view to group: ',groupName]);
	view = viewSet(view,'currentGroup',groupNum);
end

% Compute it
[tsMean,tsMedian,tsStd,tsMaxFrameDiff,tsMaxMedianDiff,tsMeanDividedByStd] = computeTimeSeriesStats(view,params);

% Make analysis structure
tsStats.name = 'timeSeriesStats';  % This can be reset by editAnalysisGUI
tsStats.type = 'timeSeriesStats';
tsStats.groupName = params.groupName;
tsStats.function = 'timeSeriesStats';
tsStats.guiFunction = 'timeSeriesStatsGUI';
tsStats.params = params;

% Install it in the view
view = viewSet(view,'newanalysis',tsStats);
view = viewSet(view,'newoverlay',tsMean);
view = viewSet(view,'newoverlay',tsMedian);
view = viewSet(view,'newoverlay',tsStd);
view = viewSet(view,'newoverlay',tsMaxFrameDiff);
view = viewSet(view,'newoverlay',tsMaxMedianDiff);
view = viewSet(view,'newoverlay',tsMeanDividedByStd);

% Save it
saveAnalysis(view,tsStats.name);

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [tsMean,tsMedian,tsStd,tsMaxFrameDiff,tsMaxMedianDiff,tsMeanDividedByStd] = ...
  computeTimeSeriesStats(view,params)

% Get nScans from view and get scanList from params
scanList = params.scanList;
nScans = viewGet(view,'nScans');

% Intialize overlay structures

% mean
tsMean.name = 'mean';
tsMean.function = 'timeSeriesStats';
tsMean.data = cell(1,nScans);
tsMean.params = params;
tsMean.colormap = jet(256);
tsMean.groupName = params.groupName;
tsMean.interrogator = 'timeSeriesStatsPlot';

% median
tsMedian.name = 'median';
tsMedian.function = 'timeSeriesStats';
tsMedian.data = cell(1,nScans);
tsMedian.params = params;
tsMedian.colormap = jet(256);
tsMedian.groupName = params.groupName;
tsMedian.interrogator = 'timeSeriesStatsPlot';

% std
tsStd.name = 'std';
tsStd.function = 'timeSeriesStats';
tsStd.data = cell(1,nScans);
tsStd.params = params;
tsStd.colormap = jet(256);
tsStd.groupName = params.groupName;
tsStd.interrogator = 'timeSeriesStatsPlot';

% max frame-to-frame diff
tsMaxFrameDiff.name = 'maxFrameDiff';
tsMaxFrameDiff.function = 'timeSeriesStats';
tsMaxFrameDiff.data = cell(1,nScans);
tsMaxFrameDiff.params = params;
tsMaxFrameDiff.colormap = jet(256);
tsMaxFrameDiff.groupName = params.groupName;
tsMaxFrameDiff.interrogator = 'timeSeriesStatsPlot';

% max diff from median
tsMaxMedianDiff.name = 'maxMedianDiff';
tsMaxMedianDiff.function = 'timeSeriesStats';
tsMaxMedianDiff.data = cell(1,nScans);
tsMaxMedianDiff.params = params;
tsMaxMedianDiff.colormap = jet(256);
tsMaxMedianDiff.groupName = params.groupName;
tsMaxMedianDiff.interrogator = 'timeSeriesStatsPlot';

% mean/std
tsMeanDividedByStd.name = 'meanDividedByStd';
tsMeanDividedByStd.function = 'timeSeriesStats';
tsMeanDividedByStd.data = cell(1,nScans);
tsMeanDividedByStd.params = params;
tsMeanDividedByStd.colormap = jet(256);
tsMeanDividedByStd.groupName = params.groupName;
tsMeanDividedByStd.interrogator = 'timeSeriesStatsPlot';

disp('Computing time series statistics...');
warning('off','MATLAB:divideByZero');
for scanIndex=1:length(scanList)
   scanNum = scanList(scanIndex);
   waitHandle = mrWaitBar(0,['Computing time-series statistics for Scan ' int2str(scanNum) ':']);

   % sliceDims: [ydim xdim] for single slice
   % volDims; [ydim xdim nslices] for single scan
   sliceDims = viewGet(view,'sliceDims',scanNum);
   volDims = viewGet(view,'dims',scanNum);

   % Initialize data with NaNs
   tsMean.data{scanNum} = NaN*ones(volDims);
   tsMedian.data{scanNum} = NaN*ones(volDims);
   tsStd.data{scanNum} = NaN*ones(volDims);
   tsMaxFrameDiff.data{scanNum} = NaN*ones(volDims);
   tsMaxMedianDiff.data{scanNum} = NaN*ones(volDims);
   tsMeanDividedByStd.data{scanNum} = NaN*ones(volDims);

   nslices = viewGet(view,'nslices',scanNum);
   for sliceNum = 1:nslices
     [tsMeanSeries,tsMedianSeries,tsStdSeries,tsMaxFrameDiffSeries,tsMaxMedianDiffSeries] = ...
       computeTimeSeriesStatsSeries(view,scanNum,sliceNum);
     tsMean.data{scanNum}(:,:,sliceNum) = reshape(tsMeanSeries,sliceDims);
     tsMedian.data{scanNum}(:,:,sliceNum) = reshape(tsMedianSeries,sliceDims);
     tsStd.data{scanNum}(:,:,sliceNum) = reshape(tsStdSeries,sliceDims);
     tsMaxFrameDiff.data{scanNum}(:,:,sliceNum) = reshape(tsMaxFrameDiffSeries,sliceDims);
     tsMaxMedianDiff.data{scanNum}(:,:,sliceNum) = reshape(tsMaxMedianDiffSeries,sliceDims);
     tsMeanDividedByStd.data{scanNum}(:,:,sliceNum) = reshape(tsMeanSeries,sliceDims)./reshape(tsStdSeries,sliceDims);
     % Update waitbar
     mrWaitBar(sliceNum/nslices,waitHandle);
   end
   mrCloseDlg(waitHandle);
end


% Fill range fields 
tsMean.range = findRange(tsMean.data);
tsMedian.range = findRange(tsMedian.data);
tsStd.range = findRange(tsStd.data);
tsMaxFrameDiff.range = findRange(tsMaxFrameDiff.data);
tsMaxMedianDiff.range = findRange(tsMaxMedianDiff.data);
tsMeanDividedByStd.range = findRange(tsMeanDividedByStd.data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [tsMeanSeries,tsMedianSeries,tsStdSeries,tsMaxFrameDiffSeries,tsMaxMedianDiffSeries] = ...
  computeTimeSeriesStatsSeries(view,scan,slice)

% Get junk frames and nframes
junkframes = viewGet(view,'junkframes',scan);
nframes = viewGet(view,'nframes',scan);

% Load tSeries
tSeries = loadTSeries(view, scan, slice);

% Reshape the tSeries
tSeries = reshapeTSeries(tSeries);

% Remove junkFrames
tSeries = tSeries(junkframes+1:junkframes+nframes,:);

tsMeanSeries = nanmean(tSeries);
tsMedianSeries = nanmedian(tSeries);
tsStdSeries = nanstd(tSeries);
tsMaxFrameDiffSeries = nanmax(abs(tSeries(2:end,:)-tSeries(1:end-1,:)));
tsMaxMedianDiffSeries = nanmax(abs(tSeries - repmat(tsMedianSeries,[nframes,1])));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function range = findRange(data)

ampMin = realmax;
ampMax = 0;
nScans = length(data);
for scan=1:nScans
  if ~isempty(data{scan})
     thisData = data{scan}(:);
     thisData = thisData(~isinf(thisData));
    ampMin = min([ampMin min(thisData)]);
    ampMax = max([ampMax max(thisData)]);
  end
end
if (ampMin <= ampMax)
  range = [ampMin ampMax];
else
  % if amp data is empty, need to make sure min < max
  range = [0 1];
end

