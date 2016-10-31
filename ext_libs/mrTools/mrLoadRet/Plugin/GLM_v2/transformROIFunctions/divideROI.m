% divideROI.m
%
%        $Id: divideROI.m 1969 2010-12-19 19:14:32Z julien $ 
%      usage: rois = divideROI(roi,numberX,numberY,numberZ,<matchNumberVoxels>)
%         by: julien besle
%       date: 27/01/2011
%
%    purpose: divide ROI into subROIs of (nearly) identical sizes 
%      input:   - numberX,numberY,numberZ" number of divisions in the X,Y,Z directions
%               - matchNumberVoxels: if 1, discard voxels to equate number of voxels per ROI (default=0)

function rois = divideROI(roi,numberX,numberY,numberZ,matchNumberVoxels)

if ~ismember(nargin,[4 5])
  help divideROI;
  return
end

if ieNotDefined('matchNumberVoxels')
  matchNumberVoxels =1;
end

if length(roi)>1
  mrWarnDlg('(divideROI) pass only one ROI at a time');
end

% %I started write a subfunction to hierarchically sort the voxels along the three dimensions,
% roi.coords = multiSort(roi.coords,[1 2 3]);
% %but the use of union in modifyROI.m ensures that the coordinates are sorted according to the first,second and last row
% and anyway there is a builtin function for that: sortrows
% Here I'll use unique, which also hierarchically sorts, because I also want to make sure there is no duplicate
roisCoords = unique(roi.coords','rows')'; 

totalNumberVoxels = size(roisCoords,2);


extraVoxels = rem(totalNumberVoxels,numberX*numberY*numberZ);
if extraVoxels
  if matchNumberVoxels
    disp(sprintf('(divideROI) Discarding %d voxels in order to divide in equal rois',extraVoxels));
    roisCoords = roisCoords(:,extraVoxels+1:end);
  else
    %randomly duplicate voxels that will be discarded later
    addedVoxels = numberX*numberY*numberZ - extraVoxels;
    roisCoords = [roisCoords roisCoords(:,ceil(totalNumberVoxels*rand(1,addedVoxels)))];
    roisCoords = sortrows(roisCoords')'; %and re-sort
  end
  totalNumberVoxels = size(roisCoords,2);
end

%Now we can start divide along the first dimension
roisCoords = reshape(roisCoords,[3 totalNumberVoxels/numberX numberX]);

%we'll sort each sub-array along the second, third  and first dimensions
for x = 1:size(roisCoords,3);
  roisCoords(:,:,x) = sortrows(roisCoords(:,:,x)',[2 3 1])'; 
end
%and divide along the second dimension
roisCoords = reshape(roisCoords,[3 totalNumberVoxels/numberX/numberY numberY numberX]);

%and same for the third dimension
for x = 1:size(roisCoords,4);
  for y = 1:size(roisCoords,3);
    roisCoords(:,:,y,x) = sortrows(roisCoords(:,:,y,x)',[3 1 2])'; 
  end
end
%divide along the third dimension and reshape into a 3D array where the third dimension is for different rois
roisCoords = reshape(roisCoords,[3 totalNumberVoxels/numberX/numberY/numberZ numberZ*numberY*numberX]);

rois = repmat(roi,1,size(roisCoords,3));
for iRoi = 1:size(roisCoords,3)
  rois(iRoi).coords = unique(roisCoords(:,:,iRoi)','rows')'; %discard duplicate voxels
  rois(iRoi).name = sprintf('%s_%03d',rois(iRoi).name,iRoi);
end



% %%%%%%%%%%5 I don't use this. 
% % it's supposed to do the same thing as sortrows
% % although I haven't finished debugging it
% function array = multiSort(array,rowOrder)
% 
% %sort along the first row
% [~,sortIndex] = sort(array(rowOrder(1),:));
% array = array(:,sortIndex);
% 
% %sort along the second and third rows
% xValues = unique(array(rowOrder(1),:));
% firstXvalue = 1;
% for i=xValues
%   %find the sub-array with identical values on the first row
%   lastXvalue = find(array(rowOrder(1),:)==i,1,'last');
%   subArray = array(:,firstXvalue:lastXvalue);
%   %sort it along the second row
%   [~,sortIndex] = sort(subArray(rowOrder(2),:));
%   subArray = subArray(:,sortIndex);
%   
%   %sort the sub array along the third rwo
%   yValues = unique(subArray(rowOrder(2),:));
%   firstYvalue = 1;
%   for j=yValues
%     %find the sub-array with identical values on the second row
%     lastYvalue = find(array(rowOrder(2),:)==j,1,'last');
%     subSubArray = subArray(:,firstYvalue:lastYvalue);
%     %sort it along the third row
%     [~,sortIndex] = sort(subSubArray(rowOrder(3),:));
%     subSubArray = subSubArray(:,sortIndex);
%     %put the sub-array back in place
%     subArray(:,firstYvalue:lastYvalue) = subSubArray;
%     firstYvalue = lastYvalue+1;
%   end
%   %put the sub-array back in place
%   array(:,firstXvalue:lastXvalue) = subArray;
%   firstXvalue = lastXvalue+1;
% end
% 
