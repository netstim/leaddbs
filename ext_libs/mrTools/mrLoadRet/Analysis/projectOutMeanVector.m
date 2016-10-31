% projectOutMeanVector.m
%
%        $Id$
%      usage: [targetROI <tSeries>] = projectOutMeanVector(v,sourceROI,targetROI,<tSeries>)
%         by: justin gardner
%       date: 05/01/08
%    purpose: computes the mean of sourceROI and then projects
%             that mean out of targetROI. Returns a structure
%             with the field tSeries which contains the
%             projected out data. If tSeries is passed
%             in then instead of loading the data from the
%             view, the data in tSeries is used. In this case, the
%             output will appear in the tSeries output.
%
%             This function works on the current scan/group in the
%             view. source/targetROI can be an roi number (or cell
%             array or roi numbers) an roi name (or cell array of
%             names), or the time series of an roi. If it is 
%             an empty matrix, all voxels in the scan will be used
%  
%             Alternatively, a params structure can be passed in which has 
%             the fields sourceROI and target ROI
%             params.sourceROI = 'left_patch';params.targetROI = [];
%             projectOutMeanVector(v,params);
%          
%             e.g. To remove the mean component calculated from l_mt from r_mt
%             r_mt = projectOutMeanVector(v,'l_mt','r_mt');
%
%             If you want to use this for two matrices without a
%             view, view can be set to []:
%
%             noise = rand(1,1000);
%             for i = 1:50
%                signal(i,:) = noise+rand(1,1000);
%             end
%             p = projectOutMeanVector([],noise,signal);
%
%            If you want to apply remove junk frames, for example when using with a scan
%            from the MotionComp group which has not already had junk frames removed, you
%            need to set optional argument 'removeJunk' to 'true':
%
%            r_mt = projectOutMeanVector(v, 'l_mt', 'r_mt', [ ],  'removeJunk', true)
% 
%    
function [targetROI tSeries] = projectOutMeanVector(v,sourceROI,targetROI,tSeries,varargin)

% check arguments
if ~any(nargin == [2 3 4 5 6])
  help projectOutMeanVector
  return
end

% Set the default (which applies to when 'projectOutMeanVector' is called
% by 'concatTSeries' or is directly applied on a scan from 'Concatenation'
% group). The default is 'not' to remove junk frames.
getArgs(varargin,{'removeJunk',false});

% This function can either be called with a params structure, or with
% two arguments specifying the source and target ROI.
params = sourceROI;
if isstruct(params)
  % get source and target ROI
  if isfield(params,'sourceROI') && isfield(params,'targetROI')
    sourceROI = params.sourceROI;
    targetROI = params.targetROI;
  else
    disp(sprintf('(projectOutMeanVector) Passed in params structure is not valid'));
    tSeries = [];
    return
  end  
else
  params = [];
  params.sourceROI = sourceROI;
  params.targetROI = targetROI;
end
if ieNotDefined('tSeries'),tSeries = [];end

% get the source and target tSeries
if isempty(params.sourceROI)
  % passed in argument is empty make an roi of the whole scan
  sourceROI = makeROIofAllVoxels(v,'source',tSeries);
elseif isnumeric(params.sourceROI) && ~isscalar(params.sourceROI)
  % passed in argument is a tSeries. Convert to an roi.
  sourceROI = makeROIfromTSeries(v,params.sourceROI,'source');
else
  % passed in argument is an roi or list of rois. So combine it
  % into a single roi and load it in
  sourceROI = loadROIList(v,params.sourceROI,tSeries);
end
if isempty(params.targetROI)
  % passed in argument is empty make an roi of the whole scan
  targetROI = makeROIofAllVoxels(v,'target',tSeries);
elseif isnumeric(params.targetROI) && ~isscalar(params.targetROI)
  % passed in argument is a tSeries. Convert to an roi.
  targetROI = makeROIfromTSeries(v,params.targetROI,'target');
else
  % passed in argument is an roi
  targetROI = loadROIList(v,params.targetROI,tSeries);
end

% get the frame numbers over which to work
% first check if this is a concat and get nFrames
if isview(v)
  concatInfo = viewGet(v,'concatInfo');
  nFrames = viewGet(v,'nFrames');
else
  concatInfo = [];
  nFrames = size(sourceROI.tSeries,2);
end
 
if isempty(concatInfo)
    
    % Check if projectOutMeanVector has been called directly on a scan from
    % MotionComp group; if so, the junk frames should be removed
    if removeJunk
        
        junkFrames = 0;
        if isview(v)
            junkFrames = viewGet(v,'junkFrames');
        end
        frameNums{1} = (junkFrames+1):(junkFrames+nFrames);
    else
        % no concatInfo, just do all frames at once
        frameNums{1} = 1:nFrames;
    end
    
else
    % get each block of frames from each scan
    % since we want to project out on a scan by scan basis
    for i = 1:concatInfo.n
        frameNums{i} = concatInfo.runTransition(i,1):concatInfo.runTransition(i,2);
    end
    disp(sprintf('(projectOutMeanVector) Scan is a concatenation. Projecting out separately for each of the %i concatenated scans',concatInfo.n));
end

% compute mean tSeries.
if sourceROI.n > 1
  meanVectorAllFrames = nanmean(sourceROI.tSeries);
else
  meanVectorAllFrames = sourceROI.tSeries;
end
% keep the ROI names
targetROI.sourceName = sourceROI.name;
% now dump sourceROI to conserve space
clear sourceROI

% cycle over each segment of a concat. If this is a single
% scan then we will be doing the whole thing in one pass
disppercent(-inf,sprintf('(projectOutMeanVector) Projecting out vector.'));
for i = 1:length(frameNums)

  % get this mean vector
  meanVector = meanVectorAllFrames(frameNums{i});
  
  % demean/detrend meanVector
  meanVector = eventRelatedDetrend(meanVector',1)';

  % make it a unit vector
  meanVector = meanVector/norm(meanVector);

  % project data on to meanVector
  % compute the magnitude of the projection of the targetROI tseries
  % on to the mean vector of the source ROI. Need to detrend/demean first
  targetROI.tSeries(:,frameNums{i}) = eventRelatedDetrend(targetROI.tSeries(:,frameNums{i})')';
  targetROI.projectionMagnitude(:,i) = targetROI.tSeries(:,frameNums{i})*meanVector';

  % now compute the normalized projection magnitude. This is the
  % projection magnitude divided by the norm of the tSeries. 
  % This computation is equivalent to the Pearson's correlation
  % coefficient which is cov(x,y)/std(x)*std(y). The reason
  % is that this in matlab format this is:
  % (x-mean(x))*(y-mean(y)'/(std(x)*std(y)*(n-1))
  % y is mean subtracted, so this gets rid of the means in the numerator
  % the x,y are divided by there norms, so the denominator
  % becomes 1. So effectively, if you take the dot product of
  % x and y which are mean subtracted and normalized: it is r!
  targetROI.r(:,i) = targetROI.projectionMagnitude(:,i)./sqrt(sum(eventRelatedDetrend(targetROI.tSeries(:,frameNums{i})',1)'.^2,2));

  % remove that component
  targetROI.tSeries(:,frameNums{i}) = targetROI.tSeries(:,frameNums{i})-targetROI.projectionMagnitude(:,i)*meanVector;

  % Calculate this magnitude with the respect to the final vector magnitude
  % so we can reconstruct original vector easily.
  targetROI.reconProjectionMagnitude(:,i) = targetROI.projectionMagnitude(:,i)./sqrt(sum(eventRelatedDetrend(targetROI.tSeries(:,frameNums{i})',1)'.^2,2));

  % if there was data passed in, then put these tSeries back into data
  if ~isempty(tSeries) && (nargout == 2)
    % special case if the roi contains all voxels
    if strcmp(targetROI.name,'allVoxels_target')
      tSeries(:,:,:,frameNums{i}) = reshape(targetROI.tSeries(:,frameNums{i}),size(tSeries,1),size(tSeries,2),size(tSeries,3),length(frameNums{i}));
    else
      % go frame by frame and copy the voxels in this roi to the output structure
      for frameNum = 1:length(frameNums{i})
	linearCoords = sub2ind(size(tSeries),targetROI.scanCoords(1,:),targetROI.scanCoords(2,:),targetROI.scanCoords(3,:),frameNums{i}(frameNum)*ones(1,targetROI.n));
	tSeries(linearCoords) = targetROI.tSeries(:,frameNum);
      end
    end
  end
  targetROI.sourceMeanVector(i,:) = meanVector;
  disppercent(i/length(frameNums));
end
disppercent(inf);

% and clear the tseries in the targetROI if we are passing back the
% data as an array
if ~isempty(tSeries) && (nargout == 2)
  targetROI.tSeries = [];
end

%% In case 'removeJunk:true' and we would like to remove the junk frames from the function outputs 
if removeJunk
    if isview(v) && isempty(concatInfo)
        
        if ~isempty(tSeries) && (nargout == 2)            
            
            tSeries = tSeries(:,:,:,frameNums{1});
        else
            targetROI.tSeries = targetROI.tSeries(:,frameNums{1});
        end
    end
end
%%
% get linear coordinates, since that is usually easier
if isview(v)
  targetROI.linearCoords = sub2ind(viewGet(v,'scanDims'),targetROI.scanCoords(1,:),targetROI.scanCoords(2,:),targetROI.scanCoords(3,:))';
end

return

mlrSmartfig('projectOutMeanVector');
% make an overlay
overlay = zeros(viewGet(v,'scanDims'));
overlay(params.targetROI{1}.linearCoords) = params.targetROI{1}.r;

% find largest 
maxMatch = first(find((max(params.targetROI{1}.r) == params.targetROI{1}.r)));
x = params.targetROI{1}.coords(1,maxMatch);
y = params.targetROI{1}.coords(2,maxMatch);
s = params.targetROI{1}.coords(3,maxMatch);
plot(tSeries(maxMatch,:));hold on
plot(params.targetROI{1}.projectionMagnitude(maxMatch)*meanVector,'r')
plot(params.targetROI{1}.tSeriesRemoved(maxMatch,:),'k');
title(sprintf('Voxel: [%i %i %i] projection mag: %0.3f',x,y,s,params.targetROI{1}.r(maxMatch)));
legend('Original','Projection component','Projection Removed');
xlabel('Time (TR)');
ylabel('MRI Signal')


keyboard



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   makeROIfromTSeries   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function roi = makeROIfromTSeries(v,tSeries,name)

roi.name = name;
roi.voxelSize = [nan nan nan];
roi.xform = eye(4);
roi.notes = sprintf('%s data passed into function projectOutMeanVector',name);
[tf roi] = isroi(roi);
roi.n = size(tSeries,1);
roi.tSeries = tSeries;
roi.coords = [];
roi.scanCoords = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   makeROIofAllVoxels   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function roi = makeROIofAllVoxels(v,name,tSeries)

% set some fields
roi.name = sprintf('allVoxels_%s',name);
roi.voxelSize = viewGet(v,'scanVoxelSize');
roi.notes = sprintf('All voxels in scan %s:%i',viewGet(v,'groupName'),viewGet(v,'curScan'));

% set the xform and xformCode to whatever is in the curent scan
sformCode = viewGet(v,'scanSformCode');
if (sformCode)
  roi.xform = viewGet(v,'scanXform');
  roi.sformCode = sformCode;
else
  roi.xform = viewGet(v,'scanQform');
  roi.sformCode = viewGet(v,'scanQformCode');
end

% add all coordinates in scan
dims = viewGet(v,'scanDims');
% the y/x flip here is intentional. it is so that the
% coordinates match when we reshape the tSeries matrix below
[y x s] = meshgrid(1:dims(2),1:dims(1),1:dims(3));
roi.coords = [reshape(x,1,prod(dims)) ; reshape(y,1,prod(dims)) ; reshape(s,1,prod(dims))];
roi.coords(4,:) = 1;

% make it into an roi
[tf roi] = isroi(roi);

roi.scanCoords = roi.coords;

% and load the data
roi.n = prod(dims);
if isempty(tSeries)
  disppercent(-inf,sprintf('(projectOutMeanVector) Loading tSeries for scan %s:%i',viewGet(v,'groupName'),viewGet(v,'curScan')));
  roi.tSeries = loadTSeries(v);
else
  disppercent(-inf,sprintf('(projectOutMeanVector) Using passed in tSeries for scan %s:%i',viewGet(v,'groupName'),viewGet(v,'curScan')));
  roi.tSeries = tSeries;
end
roi.tSeries = reshape(roi.tSeries,prod(dims(1:3)),size(roi.tSeries,4));
disppercent(inf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   thisLoadROITSeries   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function roi = thisLoadROITSeries(v,roi,tSeries)

if isempty(tSeries)
  % if no passed in time series, then just get it from the view
  roi = loadROITSeries(v,roi);
else
  % get it from the passed in tSeries
  roi = loadROITSeries(v,roi,[],[],'loadType=none');
  % get roi linear coords
  disppercent(-inf,sprintf('(projectOutMeanVector) Extracting data for roi %s from passed in data',roi.name));
  for frameNum = 1:size(tSeries,4)
    linearCoords = sub2ind(size(tSeries),roi.scanCoords(1,:),roi.scanCoords(2,:),roi.scanCoords(3,:),frameNum*ones(1,roi.n));
    roi.tSeries(:,frameNum) = tSeries(linearCoords);
    disppercent(frameNum/size(tSeries,4));
  end
  disppercent(inf);
end

%%%%%%%%%%%%%%%%%%%%%%%
%%   Load ROI List   %%
%%%%%%%%%%%%%%%%%%%%%%%
function roi = loadROIList(v,roilist,tSeries)

% make sure we are passed a cell array
roilist = cellArray(roilist);

% we will make a list of all the roi nums
roinumlist = [];

nROIs = viewGet(v,'nROIs');
for roinum = 1:length(roilist)
  % if it is a string, see if we have the roi loaded
  if isstr(roilist{roinum})
    % get the roinum
    thisroinum = viewGet(v,'roinum',roilist{roinum});
    % if it doesn't exist, then try to load it.
    if isempty(thisroinum)
      v = loadROI(v,roilist{roinum});
      thisroinum = viewGet(v,'roinum',roilist{roinum});
    end
    % now add it to the roilist if it is valid
    if ~isempty(thisroinum)
      roinumlist(end+1) = thisroinum;
    else
      disp(sprintf('(projectOutMeanVector) Could not find ROI %s. Ignoring.',roilist{roinum}));
    end
  % if it is just a roi number, make sure it is in the right range
  elseif isscalar(roilist{roinum})
    if (roilist{roinum} >=1) && (roilist{roinum} <= nROIs)
      roinumlist(end+1) = roilist{roinum};
    else
      disp(sprintf('(projectOutMeanVector) Could not find ROI %i. Ignoring.',roilist{roinum}));
    end
  end
end
 
% check to make sure we got some rois
if isempty(roinumlist)
  mrErrorDlg(sprintf('(projectOutMeanVector) No valid ROIs found'));
end

% make sure it is a unique list of roi numbers
roinumlist = unique(roinumlist);

% check to see if we have more than one roi
if length(roinumlist) > 1
  % now we will make everything into a combined ROI.
  % first get a temporary name
  tempROIName = 'Combination of';
  for roinum = 1:length(roinumlist)
    tempROIName = sprintf('%s %s',tempROIName,viewGet(v,'roiName',roinumlist(roinum)));
  end
  tempROIName = fixBadChars(tempROIName);

  v = combineROIs(v,roinumlist(1),roinumlist(1),'Union',tempROIName);

  % now, we combine the rois until we have one
  for roinum = 2:length(roinumlist)
    v = combineROIs(v,tempROIName,roinumlist(roinum),'Union');
  end
  % load the data
  roi = thisLoadROITSeries(v,tempROIName,tSeries);
  
  % and remove the temporary roi
  v = viewSet(v,'deleteroi',viewGet(v,'roinum',tempROIName));
else
  roi = thisLoadROITSeries(v,roinumlist(1),tSeries);
end  
