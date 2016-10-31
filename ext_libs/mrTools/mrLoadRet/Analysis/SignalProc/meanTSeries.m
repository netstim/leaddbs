function tseries = meanTSeries(view, groupNum, roiList, scanList, varargin)
%
% tseries = meanTseries(view, groupNum, roiList, scanList, [param1], [value1], [param2], [value2])
% 
% Computes mean tseries for each ROI and each scan, excluding junkframes
%
% groupNum: default current group
% roiList: vector of ROI numbers (default: [1:nROIs])
% scanList: vector of scan numbers (default: [1:nscans])
% other arguments as in percentTSeries  
%
% tseries: cell array (nROIs x nScans), each element of which is nFrames
% vector.
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

% Initialize result cell array
tseries = cell(length(roiList),length(scanList));

% Extract tseries for each voxel in each ROI & scan
tseriesAll = tseriesROI(view, groupNum, roiList, scanList, varargin{:});

% Loop through ROIs & scans and slices, computing means across voxels.
for iROI = 1:length(roiList)
    for iscan = 1:length(scanList)
        tseries{iROI,iscan} = nanmean(tseriesAll{iROI,iscan},2);
    end
end
return

% Test
tseries = meanTseries(MLR.views{1});

