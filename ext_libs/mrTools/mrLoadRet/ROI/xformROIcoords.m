function newcoords = xformROIcoords(coords,xform,inputVoxSize,outputVoxSize,sampRate)
%
%        $Id$
%
% newcoords = xformROIcoords(coords,xform,inputVoxSize,outputVoxSize,[sampRate])
%
% Transforms ROI coords using Xform, supersampling in each dimension to
% accumulate partial volumes, then keeping only those voxels with partial
% volumes above thresh to maintain ROI volume.
% 
% coords: 4xN matrix of coordinates (x,y,z,1).
% xform: 4x4 homogeneous transform. Note that this xform should be
% the xform from the roi to whatever coordinates you wanted. This
% function used to also add the shiftOriginXform transforms on to
% shift from 0,0,0 based to 1,1,1 but now assumes that those have
% already been composited on to the xform.
% inputVoxSize: 3-vector, size of voxels (mm) in coords
% outputVoxSize: 3-vector, size of voxels (mm) in newCoords
% sampRate: 3-vector, supersampling rate for each dimension
%           default is odd number >= 4x ratio of inputVoxSize/outputVoxSize
% 
% newcoords: 4xN matrix of (x,y,z,1) 
%
% djh, 8/98.  Modified from ROIcoords/transformROI.m in mrLoadRet-1
% 7/19/02 djh, Modified to maintain equal volumes
% 8/2005 djh, Updated to mrLoadRet-4.0

% check for empty coords
if isempty(coords)
  newcoords = [];
  return
end
% check for no transform, accounting for round off error
roundVal = 10000000;
xformRound = round(xform*roundVal)/roundVal;
inputVoxSizeRound = round(inputVoxSize*roundVal)/roundVal;
outputVoxSizeRound = round(outputVoxSize*roundVal)/roundVal;
% now check for identity xform and same voxel sizes
if isequal(xformRound,eye(4)) && isequal(inputVoxSizeRound,outputVoxSizeRound)
  % This is where coordinates get rounded - may need to change
  % this if we keep roi coordinates at finer than 1x1x1 mm resolution
  coords = round(coords);
  % get unique coordinates, do it as a linear array since it is faster
  maxCoord = repmat(max(coords(:)),1,3);
  coordsLinear = unique(mrSub2ind(maxCoord,coords(1,:),coords(2,:),coords(3,:)));
  [newcoords(1,:) newcoords(2,:) newcoords(3,:)] = ind2sub(maxCoord,coordsLinear);
  newcoords(4,:) = 1;
  return
end

if ieNotDefined('sampRate')
    if size(outputVoxSize,1) ~= 1 %if not row matrix
        outputVoxSize = outputVoxSize';
    end
  sampRate = ceil(inputVoxSize ./ outputVoxSize) .* [4,4,4];
  sampRate = 2*floor(sampRate/2) + 1;
end

if isempty(coords)
  newcoords = [];
  return;
end

% First, just transform the coords. This is insufficient because it will
% leave holes if the input voxels are larger than the output voxels.
%
tmpNewCoords = xform * coords;

% Find bounding (min and max) volume.
%
minPos = [min(tmpNewCoords(1,:));min(tmpNewCoords(2,:));min(tmpNewCoords(3,:));1];
maxPos = [max(tmpNewCoords(1,:));max(tmpNewCoords(2,:));max(tmpNewCoords(3,:));1];
minPos = floor(minPos)-[round((inputVoxSize ./ outputVoxSize)+1) 0]'; %JB: this is to make sure that all potentially transformed coordinates fall in the bounding box
maxPos = ceil(maxPos)+[round((inputVoxSize ./ outputVoxSize)+1) 0]';  % it used to be +/-[2 2 2 0], but this doesn't work if new voxel size < old_voxel_size/2
dims = (maxPos(1:3)-minPos(1:3)+ones(3,1))';

% Initialize accumulator for partial volume calculation, a vector
% of length appropriate to index the bounding volume.
%
accum = zeros(1,prod(dims));

% Calculate offsets that will be added within the loop to do the
% partial voluming.
%
xoffsets=[-.5+1/(2*sampRate(1)):1/sampRate(1):.5-1/(2*sampRate(1))];
yoffsets=[-.5+1/(2*sampRate(2)):1/sampRate(2):.5-1/(2*sampRate(2))];
zoffsets=[-.5+1/(2*sampRate(3)):1/sampRate(3):.5-1/(2*sampRate(3))];
% xoffsets=[0:1/sampRate(1):1-1/sampRate(1)];
% yoffsets=[0:1/sampRate(2):1-1/sampRate(2)];
% zoffsets=[0:1/sampRate(3):1-1/sampRate(3)];

% Divide alpha by prod(sampRate) to get partial volume for the
% supersampled voxels.
%
alpha = repmat(1/prod(sampRate),[1 size(coords,2)]);

% Loop through supersamples, transform them, and accumulate
% partial volume.
%
for ioff=1:length(xoffsets)
  xoff=xoffsets(ioff);
  for yoff=yoffsets
    for zoff=zoffsets
      % Add offset
      tmpNewCoords(1:3,:) = coords(1:3,:) + repmat([xoff;yoff;zoff],[1,size(coords,2)]);
      % Transform
      tmpNewCoords = xform * tmpNewCoords;
      % Round and subtract minPos
      tmpNewCoords(1:3,:) = round(tmpNewCoords(1:3,:)) - repmat(minPos(1:3),[1,size(tmpNewCoords,2)]);
%          % jg: make sure that tmpNewCoords doesn't go to zero--this      % JB  i don't think we need this here anymore (see above)
%          % seems to happen because of a rounding error from the above    %
%          % statement.                                                    %
%          if sum(tmpNewCoords(:)==0)                                      %   
%       %disp(sprintf('(xformROIcoords) Zero index corrected'));           %
%       tmpNewCoords((tmpNewCoords(:)==0)) = 1;                            %
%          end                                                             %
       % Convert to indices
      indices = sub2ind(dims,tmpNewCoords(1,:),tmpNewCoords(2,:),tmpNewCoords(3,:));%JB: not sure why that would be necessary... I thought sub2ind returned an error if coordinates outside dims
      indices = indices(~isnan(indices));
      % Accumulate partial volume. Need to do it in a loop
      % instead of:
      %    accum(indices) = accum(indices) + alpha;
      % because an index can appear twice in indices and we want
      % to accumulate them both.
      for jj=1:length(indices)
	accum(indices(jj)) = accum(indices(jj)) + alpha(jj);
      end
    end
  end
end

% Build newROIcoords
%
%voxels in new space are sorted according to how many old-space supersampled voxels they received
[sortedAccum,indices] = sort(accum);
nonZeroSize = nnz(accum);
%Compute number of voxels necessary to maintain the ROI volume with the new voxel size
newROIsize = round(prod(inputVoxSize)*size(coords,2) / prod(outputVoxSize));
%we don't want to keep voxels that have zero partial voluming
newROIsize = min(nonZeroSize,newROIsize);
%the volume is conserved by keeping the appropriate number of voxels with highest partial voluming
indices = indices(length(indices)-newROIsize+1:length(indices));
if ~isempty(indices)
  [newcoords1,newcoords2,newcoords3] = ind2sub(dims,indices);
  newcoords = [newcoords1; newcoords2; newcoords3; ones(size(newcoords1))];
  newcoords(1:3,:) = newcoords(1:3,:) + repmat(minPos(1:3),[1,length(indices)]);
else
  newcoords = [];
end

% this is just here to check the coordinates created by
% not going through this logic, but just returing what
% was there
%if ~ieNotDefined('testnewcoords')
%  a = sort(mrSub2ind([256 256 256],newcoords(1,:),newcoords(2,:),newcoords(3,:)));
%  b = sort(mrSub2ind([256 256 256],testnewcoords(1,:),testnewcoords(2,:),testnewcoords(3,:)));
%  disp(sprintf('Consistency check: %i',isequal(a,b)));
%end

return;

%%%%%%%%%%%%%%
% Debug/test %
%%%%%%%%%%%%%%

coords1 = [0; 0; 0; 1];
coords2 = [1 2 3 4;
	   1 1 1 1;
	   1 1 1 1;
	   1 1 1 1];
xform = [1 0 0 0;
	 0 1 0 0;
	 0 0 1 0.5;
	 0 0 0 1];
inputVoxSize = [1,1,1];
outputVoxSize = [1,1,1/2];
xformROIcoords(coords1,xform,inputVoxSize,outputVoxSize)
xformROIcoords(coords2,xform,inputVoxSize,outputVoxSize)
