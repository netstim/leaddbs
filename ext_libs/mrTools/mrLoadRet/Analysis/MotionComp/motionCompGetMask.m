% motionCompGetmask.m
%
%      usage: motionCompGetmask()
%         by: justin gardner
%       date: 03/02/10
%    purpose: Gets a mask to mask out the tSeries with from the
%             params returned by motionCompGUImrParams
%
function [mask view] = motionCompGetMask(view,params,scanNum,groupNum)

maskROINum = [];
% get the maskROINum if we are using an ROI to mask the tSeries
if isfield(params,'useMask') && params.useMask
  % get the roi num
  maskROINum = viewGet(view,'roiNum',params.maskROI);
  % if roi is not already loaded, then load it and get num
  if isempty(maskROINum)
    view = loadROI(view,params.maskROI);
    maskROINum = viewGet(view,'roiNum',params.maskROI);
  end
  % couldn't find roi, print warning and continue
  if isempty(maskROINum),disp(sprintf('(motionComp) Could not load mask %s',params.maskROI));end
end

% get scan nums
scanDims = viewGet(view,'scanDims',scanNum,groupNum);

% get the mask to use, if any
if ~isempty(maskROINum)
  % if we have an roi, then we only use voxels within the mask ROI
  maskCoords = getROICoordinates(view,maskROINum,scanNum,groupNum);
  maskLinearCoords = sub2ind(scanDims,maskCoords(1,:),maskCoords(2,:),maskCoords(3,:));
  mask = zeros(scanDims);
  mask(maskLinearCoords) = 1;
else
  % otherwise, we use all values (i.e. like no mask)
  mask = ones(scanDims);
end

