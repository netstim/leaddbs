% convertROICorticalDepth.m
%
%        $Id$
%      usage: convertROICorticalDepth()
%         by: justin gardner
%       date: 10/15/07
%    purpose: used to extend or restrict ROI coordinates across
%             cortical depths
%
%  12/8/08 Modified by Taosheng Liu to take params. If params is set, GUI
%  will not show for setting params, also it assumes then all ROIs
%  associated with a view will be converted.

function [v params] = convertROICorticalDepth(v,params,varargin)

% check arguments
if ~any(nargin == [1 2 3 4])
  help convertROICorticalDepth
  return
end

eval(evalargs(varargin,[],[],{'justGetParams','defaultParams'}));
if ieNotDefined('justGetParams'),justGetParams = 0;end
if ieNotDefined('defaultParams'),defaultParams = 0;end
% number of rois
numrois = viewGet(v,'numberofrois');
if numrois == 0
  mrWarnDlg('(convertROICorticalDepth) No currently loaded ROIs');
  return
end
roinames = viewGet(v,'roiNames');

if ieNotDefined('params')
  % get cortical depth
  corticalDepth = viewGet(v,'corticalDepth');
  referenceDepth= mean(corticalDepth);
  if length(corticalDepth)==2 && corticalDepth(1)~=corticalDepth(2)
    minDepth = corticalDepth(1);
    maxDepth = corticalDepth(2);
  else
    minDepth = 0;
    maxDepth = 1;
  end
  depthStep = 1/(viewGet(v,'corticalDepthBins')-1);
  incdecString = sprintf('incdec=[-%f %f]',depthStep,depthStep);
  paramsInfo = {};
  paramsInfo{end+1} = {'conversionType',{'Project through depth','Restrict to reference depth'},'type=popupmenu','If you set project through depth, then this will add all the voxels from each cortical depth that are in the same position as the ones at the reference depth. If you set to restrict to reference depth, this will remove any voxels that are not on the reference depth (note that you will still see some voxels on other depths, but those are voxels that exist at the reference depth--also, voxels that do not exist on this flat map will not be affected)'};
  paramsInfo{end+1} = {'referenceDepth',referenceDepth,'min=0','max=1',incdecString,'The cortical depth to start from'};
  paramsInfo{end+1} = {'minDepth',minDepth,'min=0','max=1',incdecString,'The start depth'};
  paramsInfo{end+1} = {'depthStep',depthStep,'min=0','max=1',incdecString,'The depth step (i.e. we will go from minDepth:depthStep:maxDepth (skipping the reference depth), including or excluding voxels'};
  paramsInfo{end+1} = {'maxDepth',maxDepth,'min=0','max=1',incdecString,'The end depth'};
  paramsInfo{end+1} = {'excludeOtherVoxels',1,'type=checkbox','If ROI voxels exist oustide the projected surface, they will be remove. Uncheck to keep them. this option is ignored if restriction is selected'};
  if defaultParams
    params = mrParamsDefault(paramsInfo);
  else
    % put up some parameter choices
    params = mrParamsDialog(paramsInfo,'ROI cortical depth conversion');
  end
  % now select rois
  % put up a dialog with rois to select
  paramsDialog = {};
  for roinum = 1:length(roinames)
    helpinfo = sprintf('Convert cortical depth of ROI %i: %s',roinum,roinames{roinum});
    paramsDialog{end+1} = {fixBadChars(roinames{roinum}),0,'type=checkbox',helpinfo};
  end
  paramsDialog{end+1} = {'all',0,'type=checkbox','Select all ROIs'};
  % put up dialog
  whichROI = mrParamsDialog(paramsDialog,sprintf('Select ROIs to convert cortical depth'));
else
  disp('(convertROICorticalDepth) coverting all ROIs in the view');
  whichROI.all=1;
end
if isempty(params),return,end
% just return parameters
if justGetParams, return, end

currentROI = viewGet(v,'currentROI');
% now go through and do conversion
if ~isempty(whichROI)
  needToRefresh = 0;
  % now go through and convert anything the user selected
  for roinum = 1:length(roinames)
    if whichROI.all || whichROI.(fixBadChars(roinames{roinum})) 
      needToRefresh = 1;
      disppercent(-inf,sprintf('(convertROICorticalDepth) Processing ROI %i:%s',roinum,roinames{roinum}));
      % get the roi
      v = viewSet(v,'curROI',roinum);
      % now try to figure out what base this was created on
      roiCreatedOnBase = viewGet(v,'roiCreatedOnBase',roinames{roinum});
      if isempty(roiCreatedOnBase)
	disp(sprintf('(convertROICorticalDepth) Converting %s based on base:%s because roiCreatedOnBase has not been set.',roinames{roinum},viewGet(v,'baseName')));
	baseNum = viewGet(v,'curBase');
      else
	% get the basenumber for the base that this was created on
	baseNum = viewGet(v,'baseNum',roiCreatedOnBase);
	if isempty(baseNum)
	  disp(sprintf('(convertROICorticalDepth) Converting %s based on base:%s because base:%s which this roi was created on is not loaded',roinames{roinum},viewGet(v,'baseName'),roiCreatedOnBase));
	  baseNum = viewGet(v,'curBase');
	end
      end
      % get the roi transformation in order to set the coordinates later
      base2roi = viewGet(v,'base2roi',roinum,baseNum);
      % get the roiBaseCoords
      roiBaseCoords = getROICoordinates(v,roinum,[],[],'baseNum',baseNum);
      if isempty(roiBaseCoords)
        disppercent(inf);
        mrWarnDlg(sprintf('(convertROICorticalDepth) %s has no coordinates on this flat',roinames{roinum}));
        continue;
      end
      % get base info
      baseVoxelSize = viewGet(v,'baseVoxelSize',baseNum);
      baseCoordMap = viewGet(v,'baseCoordMap',baseNum,params.referenceDepth);
      baseDims = baseCoordMap.dims;
      baseCoordMap = round(baseCoordMap.coords);
      referenceBaseCoordMap = mrSub2ind(baseDims,baseCoordMap(:,:,:,1),baseCoordMap(:,:,:,2),baseCoordMap(:,:,:,3));
      referenceBaseCoordMap = referenceBaseCoordMap(:);
      % get roi linear coordinates
      roiBaseCoordsLinear = mrSub2ind(baseDims,roiBaseCoords(1,:),roiBaseCoords(2,:),roiBaseCoords(3,:));
      % now find which baseCoords are in the current roi
      [isInROI roiInBase] = ismember(referenceBaseCoordMap,roiBaseCoordsLinear);
      % get the roi base coordinates that are found in base
      roiInBase = unique(setdiff(roiInBase,0));
      % if we don't find most of the coordinates, then
      % probably good to complain and give up
      if (length(roiInBase)/length(roiBaseCoordsLinear)) < 0.1
        disppercent(inf);
        mrWarnDlg(sprintf('(convertROICorticalDepth) !!! %s has less than %0.0f%% coordinates on surface %s. Perhaps you need to load the base that it was orignally created on. !!!',roinames{roinum},ceil(100*(length(roiInBase)/length(roiBaseCoordsLinear))),viewGet(v,'baseName',baseNum)));
        continue;
      end
      % make sure to keep the voxels at the reference depth
      roiBaseCoordsReferenceLinear = roiBaseCoordsLinear(ismember(roiBaseCoordsLinear,referenceBaseCoordMap));
      
      if strcmp(params.conversionType,'Project through depth')
        %clear all voxels if we're not keeping voxels outside the projection
        if params.excludeOtherVoxels
          % remove everything from the ROI
          roiBaseCoords(4,:) = 1;
          v = modifyROI(v,roiBaseCoords,base2roi,baseVoxelSize,0);
        end
        roiBaseCoordsLinear=[];
        % now get each cortical depth, and add/remove voxels
        corticalDepths = params.minDepth:params.depthStep:params.maxDepth;
        baseCoordMap = viewGet(v,'baseCoordMap',baseNum,corticalDepths);
        for iDepth = 1:size(baseCoordMap.coords,5)
          % get the coordinates at this depth
          baseCoords = round(baseCoordMap.coords(:,:,:,:,iDepth));
          baseCoords = mrSub2ind(baseDims,baseCoords(:,:,:,1),baseCoords(:,:,:,2),baseCoords(:,:,:,3));
          baseCoords = baseCoords(:);
          % add the coordinates to our list
          roiBaseCoordsLinear = union(roiBaseCoordsLinear,baseCoords(isInROI));
        end
        roiBaseCoordsLinear = roiBaseCoordsLinear(~isnan(roiBaseCoordsLinear));
        % now convert back to regular coords
        roiBaseCoords = [];
        [roiBaseCoords(1,:) roiBaseCoords(2,:) roiBaseCoords(3,:)] = ind2sub(baseDims,roiBaseCoordsLinear);
        roiBaseCoords(4,:) = 1;
        % add them to the ROI
        v = modifyROI(v,roiBaseCoords,base2roi,baseVoxelSize,1);
      else
        % get current coords
        curROICoords = viewGet(v,'roiCoords',roinum);
        % remove them from the ROI
        v = modifyROI(v,roiBaseCoords,base2roi,baseVoxelSize,0);
        % but make sure we have the voxels at the reference depth
        roiBaseCoords = [];
        [roiBaseCoords(1,:) roiBaseCoords(2,:) roiBaseCoords(3,:)] = ind2sub(baseDims,roiBaseCoordsReferenceLinear);
        roiBaseCoords(4,:) = 1;
        v = modifyROI(v,roiBaseCoords,base2roi,baseVoxelSize,1);
        % and save for undo (note we do this instead of allowing
        % modifyROI to do it since we have called modifyROI twice)
        v = viewSet(v,'prevROIcoords',curROICoords);
      end
      disppercent(inf);
    end
  end
  v = viewSet(v,'currentROI',currentROI);
  if needToRefresh
    refreshMLRDisplay(viewGet(v,'viewNum'));
  end
end

