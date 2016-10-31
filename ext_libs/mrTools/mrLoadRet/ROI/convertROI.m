
% convertROI.m
%
%        $Id$
%      usage: convertROI()
%         by: justin gardner
%       date: 10/15/07
%    purpose: convert ROIs to current base anatomy/scan xform/voxelSize
%
function thisView = convertROI(thisView)

% check arguments
if ~any(nargin == [1])
  help convertROI
  return
end

%get number of loaded bases
numBases = viewGet(thisView,'numberofBaseVolumes');
curBaseNum = viewGet(thisView,'currentBase');
for iBase = 1:numBases
% base xform and voxel size
   baseXform{iBase} = viewGet(thisView,'baseXform',iBase);
%   baseVoxelSize{iBase} = viewGet(thisView,'baseVoxelSize',iBase);
   baseName{iBase} = viewGet(thisView,'baseName',iBase);
end
curScanXform = viewGet(thisView,'scanxform');

% number of rois
numRois = viewGet(thisView,'numberofrois');
% put up a dialog with rois to delete
roiNames = viewGet(thisView,'roiNames');
defaultSpaceName = 'Unknown/Not loaded';
roiSpace = mat2cell(repmat(defaultSpaceName,numRois,1),ones(numRois,1),length(defaultSpaceName));
paramsDialog = {};
paramsDialog{end+1} = {'conversionType',{'Convert coordinates','Only adopt xform'},'type=popupmenu',...
   'Convert will convert the ROI coordinates to the destination space. This is not normally necessary (as ROIs are always converted on the fly to the base anatomy), but you may want to do this if you originally defined the ROI on a low resolution scan and want to have the ROI represented in a higher resolution or vice-vers. ''Adopt xform'' is only necessary in even more rare cases. This does not convert the roi voxels, but simply adopts the xform and voxel size of the base anatomy/scan. If for example you changed the qform/sform of your base anatomy/scan you might need to use this.'};
paramsDialog{end+1} = {'destinationSpace',[baseName,{'current scan'}],'type=popupmenu',...
   'which coordinate space you want to convert the ROI to/adopt the xform from'};
for roinum = 1:numRois
  
  % roi name
  helpinfo = sprintf('Convert ROI %i: %s',roinum,roiNames{roinum});
  paramsDialog{end+1} = {fixBadChars(roiNames{roinum}),0,'type=checkbox',helpinfo};
  
  % roi xform and voxel size
  roiXform = viewGet(thisView,'roiXform',roinum);
  if isequal(roiXform,curScanXform)
     roiSpace{roinum} = 'Current Scan';
  elseif isequal(roiXform,baseXform{curBaseNum})
     roiSpace{roinum} = ['Current Base (' baseName{curBaseNum} ')'];
  else         %here only find first base with identical xform, could be modified to find all bases
    for iBase = 1:numBases
      if isequal(baseXform{iBase},roiXform)
         roiSpace{roinum} = baseName{iBase};
         break;
      end
    end
    %if we're still here, that means we haven't found the space
    %display the voxel size to give a clue to user
    roiVoxelSize = viewGet(thisView,'roiVoxelSize',roinum);
    paramsDialog{end+1} = {sprintf('%s_voxelSize',fixBadChars(roiNames{roinum})),roiVoxelSize,'editable=0',sprintf('Current voxel size for roi %s',roiNames{roinum})};
  end
  paramsDialog{end+1} = {[fixBadChars(roiNames{roinum}) ' space'],roiSpace{roinum},'editable=0',['Current coordinate space for roi ' roiNames{roinum}]};
end

% put up dialog
  %params = mrParamsDialog(paramsDialog,sprintf('Select ROIs to convert to [%0.2g %0.2g %0.2g] resolution',baseVoxelSize(1),baseVoxelSize(2),baseVoxelSize(3)));
  params = mrParamsDialog(paramsDialog,'Select ROIs to convert');

if isempty(params) %user pressed cancel
   return;
end



if strcmp(params.destinationSpace,'current scan')
      newXform = curScanXform;
      newSformCode = viewGet(thisView,'scanSformCode');
      newVol2mag = viewGet(thisView,'scanVol2mag');
      newVol2tal = viewGet(thisView,'scanVol2tal');
      newVoxelSize = viewGet(thisView,'scanVoxelSize');
      whichVolume = [];
      baseNum = [];
      
else
      baseNum = find(ismember(baseName,params.destinationSpace));
      newXform = baseXform{baseNum};
      newSformCode = viewGet(thisView,'baseSformCode',baseNum);
      newVol2mag = viewGet(thisView,'baseVol2mag',baseNum);
      newVol2tal = viewGet(thisView,'baseVol2tal',baseNum);
      newVoxelSize = viewGet(thisView,'baseVoxelSize',baseNum);
      whichVolume = 0;
      
end

currentROIName = viewGet(thisView,'roiname');
needToRefresh = 0;
% now go through and do conversion
for roinum = 1:length(roiNames)
 if isfield(params,fixBadChars(roiNames{roinum})) && params.(fixBadChars(roiNames{roinum}))
   % get the roi
   thisroinum = viewGet(thisView,'roinum',roiNames{roinum});
   roi = viewGet(thisView,'ROI',thisroinum);
   if strcmp(params.conversionType,'Convert coordinates')
     disp(sprintf('(convertROI) Converting ROI %i:%s to %s coordinate space',roinum,roiNames{roinum},params.destinationSpace));
     if ~isempty(roi.coords)
         roiBaseCoords = getROICoordinates(thisView,thisroinum,whichVolume,[],'baseNum',baseNum);
         roi.coords = roiBaseCoords;
      else
         mrWarnDlg(sprintf('(convertROI) ROI %i:%s has empty coordinates in transformation, skipping conversion...',roinum,roiNames{roinum}));
      end
   else
     disp(sprintf('(convertROI) Adopting %s xform for ROI %i:%s',params.destinationSpace,roinum,roiNames{roinum}));
   end
   % set roi xform and voxel size
   roi.sformCode = newSformCode;
   roi.xform = newXform;
   % now set the vol2mag and vol2tal correctly
   roi.vol2mag = newVol2mag;
   roi.vol2tal = newVol2tal;
   roi.voxelSize = newVoxelSize;
   % remove the roi
   %thisView = viewSet(thisView,'deleteROI',thisroinum);   %JB: do not remove the roi. this way, user is prompted to overwrite the file of change the name
   % and add it back with the new coordinates
   thisView = viewSet(thisView,'newROI',roi);
   needToRefresh = 1;
 end
end
if needToRefresh
 % select the same current roi
 thisView = viewSet(thisView,'curROI',viewGet(thisView,'roinum',currentROIName));
 refreshMLRDisplay(viewGet(thisView,'viewNum'));
end

