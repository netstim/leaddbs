% mrReduceVolume.m
%
%        $Id$
%      usage: mrReduceVolume()
%         by: julien besle
%       date: 10/2010
%    purpose: reduces volume to a minimum voxel size
%

function [reducedVolume, reductionMatrix] = mrReduceVolume(volume,voxelSize,minVoxelSize, imageName)

if ieNotDefined('imageName')
   imageName = '';
end

% reduce volume if voxelSize <minVoxelSize to reduce memory requirements
if any(voxelSize<minVoxelSize)
   reductionFactor = ceil(max(1,minVoxelSize./voxelSize));
   wbh = mrMsgBox([ 'Please wait: ' imageName ' voxel size is ' mat2str(voxelSize) '. Reducing to ' mat2str(voxelSize.*reductionFactor) '...']);
   lpf=[.0625 .25 .375 .25 .0625]';
   inpFilt = convXYZsep(volume, lpf, lpf, lpf, 'repeat', 'same');
   reducedVolume = inpFilt(1:reductionFactor(1):size(inpFilt,1),1:reductionFactor(2):size(inpFilt,2),1:reductionFactor(3):size(inpFilt,3));
   %reductionMatrix is the matrix that transforms the coordinates of the reduced volume to the coordinates of the original volume
   reductionMatrix = diag([reductionFactor; 1]);
   mrCloseDlg(wbh);
else
   reducedVolume = volume;
   reductionMatrix = eye(4);
end


