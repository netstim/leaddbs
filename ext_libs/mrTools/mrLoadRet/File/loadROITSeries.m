function rois = loadROITSeries(view,roiname,scanList,groupNum,varargin)
% loadROITSeries.m
%
%      usage: rois = loadROITSeries(view,<roiname>,<scanList>,<groupNum>,<varargin>)
%         by: justin gardner
%       date: 03/22/07
%    purpose: load the time series for a roi, without roiname
%             specified, brings up selection dialog. roiname
%             may be a cell array. scanList and groupNum
%             default to current scan/group. If roiname is a roi
%             struct instead of a name then it will use that roi
%             instead of loading the roi from disk. Also, if roiname
%             is a number or cell array of numbers then it will use
%             the corresponding ROI from the view. 
%        e.g.:
%
% v = newView
% rois = loadROITSeries(v,[],1,1);
%
%             to load the roi coordinates, but not the time series
%
% v = newView
% rois = loadROITSeries(v,[],1,1,'loadType=none');
%
%             If a time series has nan points in it, this function
%             will drop that voxels time series (this usually happens
%             because the ROI contains a voxel that is on the edge of
%             the scan and the motion correction has brought in a nan
%             value). If you want to return all time series, regardless
%             of whether there is a nan voxel or not, pass in the
%             optional argument, keepNAN=true
%
% rois = loadROITSeries(v,[],1,1,'keepNAN',true);
%
%             You can also get voxels for the ROI that match the voxels
%             from a separate scan (sometimes useful for classification. See
%             the function getROICoordinatesMatching for more info)
%
% rois = loadROITSeries(v,[],1,1,'keepNAN',true,'matchScanNum=1','matchGroupNum=2');
%
%             You can also use the "straightXform" (without doing the fancy matching
%             of ROI size) see getROICoordinates for more info. Defaults to 0.
%
% rois = loadROITSeries(v,[],1,1,'keepNAN',true,'straightXform=1');
%   
%
% see also tseriesROI
rois = {};

% check arguments
if nargin < 1
    help loadROITSeries
    return
end

% no view specified
if ieNotDefined('view')
    view = newView;
end
if ~isview(view)
    disp(sprintf('(loadROITSeries) Passed in view is not valid'));
    return
end
% get the roi directory
roidir = viewGet(view,'roidir');

% get group and scan
if ieNotDefined('groupNum')
  groupNum = viewGet(view,'currentGroup');
end
groupName = viewGet(view,'groupName',groupNum);
if ieNotDefined('scanList')
  scanList = viewGet(view,'currentScan');
end

% set the current group
view = viewSet(view,'currentGroup',groupNum);

% if there is no roi, ask the user to select
if ieNotDefined('roiname')
  roiname = mlrGetPathStrDialog(viewGet(view,'roiDir'),'Choose one or more ROIs','*.mat','on');
end

%make into a cell array
roiname = cellArray(roiname);

% set the default arguments. loadType
% sets the way to load the time series
% possible values are 'vox' which loads each
% time series voxel by voxel (slow but less memory
% intensive), or 'block' which loads the block
% of the image around the ROI and then subselects
% the voxels needed (default--fast). Note, both
% block and vox will return the same voxel time
% series, they just differ in how they access the data
% from disk. Set to 'none' to not load the time series.
getArgs(varargin,{'loadType=block','keepNAN',false,'matchScanNum=[]','matchGroupNum=[]','straightXform=0'});

% if user has asked for a match roi (that is, the ROIs should be created with
% coordinates that are a voxel for voxel match with the passed in scan number.
% this is useful for classification protocols. Instead of calling getROICoordinates
% this function will use getROICoordinatesMatching instead.
if ~isempty(matchScanNum) && isempty(matchGroupNum)
  matchgroupNum = groupNum;
end


% load the rois in turn
for roinum = 1:length(roiname)
    % if the name is a string which is an already loaded roi
    if ischar(roiname{roinum}) & ~isempty(viewGet(view,'roinum',getLastDir(roiname{roinum})))
        roiname{roinum} = viewGet(view,'roinum',roiname{roinum});
    end
    % see if we have to paste roi directory on
    if ischar(roiname{roinum}) && ~isfile(sprintf('%s.mat',stripext(roiname{roinum})))
        roiname{roinum} = fullfile(roidir,stripext(roiname{roinum}));
    end
    % check for file
    if ischar(roiname{roinum}) && ~isfile(sprintf('%s.mat',stripext(roiname{roinum})))
        disp(sprintf('(loadROITSeries) Could not find roi %s',roiname{roinum}));
        dir(fullfile(roidir,'*.mat'))
    elseif isnumeric(roiname{roinum}) && ((roiname{roinum} < 1) || (roiname{roinum} > viewGet(view,'numberOfROIs')))
        disp(sprintf('(loadROITSeries) No ROI number %i (number of ROIs = %i)',roiname{roinum},viewGet(view,'numberOfROIs')));
    else
        % load the roi, if the name is actually a struct
        % then assume it is an roi struct. if it is a number choose
        % from a loaded roi
        if ischar(roiname{roinum})
            roi = load(roiname{roinum});
        elseif isnumeric(roiname{roinum})
            thisroi = viewGet(view,'roi',roiname{roinum});
            clear roi;
            roi.(fixBadChars(thisroi.name)) = thisroi;
        else
            roi.(fixBadChars(roiname{roinum}.name)) = roiname{roinum};
        end
        roiFieldnames = fieldnames(roi);
        % get all the rois
        for roinum = 1:length(roiFieldnames)
            for scanNum = 1:length(scanList)
                % get current scan number
                scanNum = scanList(scanNum);
                rois{end+1} = roi.(roiFieldnames{roinum});
                % set a field in the roi for which scan we are collecting from
                rois{end}.scanNum = scanNum;
                rois{end}.groupNum = groupNum;
                % convert to scan coordinates
		if isempty(matchScanNum)
		  rois{end}.scanCoords = getROICoordinates(view,rois{end},scanNum,groupNum,'straightXform',straightXform);
		else
		  rois{end}.scanCoords = getROICoordinatesMatching(view,rois{end},scanNum,matchScanNum,groupNum,matchGroupNum);
		end
		  
                % if there are no scanCoords then set to empty and continue
                if isempty(rois{end}.scanCoords)
                    rois{end}.n = 0;
                    rois{end}.tSeries = [];
                    continue;
                end
                % get x y and s in array form
                x = rois{end}.scanCoords(1,:);
                y = rois{end}.scanCoords(2,:);
                s = rois{end}.scanCoords(3,:);
                % set the n
                rois{end}.n = length(x);
                % load the tseries, voxel-by-voxel
                disppercent(-inf,sprintf('(loadROITSeries) Loading tSeries for %s from %s: %i',rois{end}.name,groupName,scanNum));
		% for now we always load by block, but if memory is an issue, we can
                % switch this if statement and load voxels indiviudally from file
                if strcmp(loadType,'vox')
                    % load each voxel time series indiviudally
                    for voxnum = 1:rois{end}.n
                        rois{end}.tSeries(voxnum,:) = squeeze(loadTSeries(view,scanNum,s(voxnum),[],x(voxnum),y(voxnum)));
                        disppercent(voxnum/rois{end}.n);
                    end
                    disppercent(inf);
                elseif strcmp(loadType,'block');
                    % load the whole time series as a block (i.e. a block including the min and max voxels)
                    % this is usually faster then going back and loading each time series individually
                    % but is more inefficient with memory
                    tSeriesBlock = loadTSeries(view,scanNum,[min(s) max(s)],[],[min(x) max(x)],[min(y) max(y)]);
		    % preallocate memory for tSeries
		    rois{end}.tSeries = nan(rois{end}.n,size(tSeriesBlock,4));
                    % now go through and pick out the voxels that we need.
                    for voxnum = 1:rois{end}.n
                        rois{end}.tSeries(voxnum,:) = squeeze(tSeriesBlock(x(voxnum)-min(x)+1,y(voxnum)-min(y)+1,s(voxnum)-min(s)+1,:));
                        disppercent(voxnum/rois{end}.n);
                    end
                    clear tSeriesBlock;
                    disppercent(inf);
                else
                    disppercent(inf);
                    disp(sprintf('(loadROITSeries) Not loading time series (loadType=%s)',loadType));
                end
                if ~strcmp(loadType, 'none')
                  nanVoxels=sum(isnan(rois{end}.tSeries),2)>0;  %returns 1 for a voxel if any number in its timeseries is nan
                  if sum(nanVoxels>0)
		    if keepNAN
		      % just display a message indicating how many nan voxels we have
		      disp(sprintf('(loadROITSeries) %s %s:%i has %i/%i (%0.2f%%) voxels with NaNs in timeseries',rois{end}.name,viewGet(view,'groupName',groupNum),scanNum,sum(nanVoxels),rois{end}.n,100*sum(nanVoxels)/rois{end}.n));
		    else
		      % remove the nan voxels from the roi and report what we did
		      rois{end}.tSeries=rois{end}.tSeries(~nanVoxels,:);
		      rois{end}.scanCoords=rois{end}.scanCoords(:,~nanVoxels');
		      disp(sprintf('(loadROITSeries) Deleted %i/%i (%0.2f%%) voxels with NaNs in timeseries in ROI %s %s:%i.',sum(nanVoxels),rois{end}.n,100*sum(nanVoxels)/rois{end}.n,rois{end}.name,viewGet(view,'groupName',groupNum),scanNum));
		      rois{end}.n=sum(~nanVoxels);
		    end
                  end
                end
            end
        end
    end
end
if length(rois) == 1
    rois = rois{1};
end

