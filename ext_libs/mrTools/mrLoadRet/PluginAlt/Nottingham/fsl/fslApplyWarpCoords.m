% warpedCoords = fslApplyWarpCoords(coords, coordsVoxelSize, warpResolution, warpCoefFilename, tempFilename, hdr, verbose)
%
%   applies non-linear FSL registration to coordinates
%
%  input: 
%          coords: coordinates in the space of the input FNIRT volume, in homogeneous format (4*n matrix, with 1s on the last row)
% coordsVoxelSize: voxel size of the space of the coordinates
%  warpResolution: resolution in mm of the warp field computed using fnirtfileutils. don't make too small because uses too much memory
%                   better to make it .5 or 1 mm, because more efficient to interpolate using interp3 than computing high-resolution warp field volume
%    tempFilename: name of the temporary file to write the data to the disc
%             hdr: header of a volume in the space of the input FNIRT volume for dimensions and voxel size (transformation matrix is not taken into account, (but qform might determine the voxel size... ?)
%
% jb 16/04/2011
%
% $Id: fslApplyWarp.m 2107 2011-04-17 19:49:52Z julien $

function warpedCoords = fslApplyWarpCoords(coords,coordsVoxelSize,warpResolution, inversewarpfile, hdr)


if ieNotDefined('verbose')
   verbose = true;
end

 scalingFactor = round(hdr.pixdim(2:4)'./warpResolution);
    dataSize = hdr.dim(2:4)'.*scalingFactor;
    data = ones(dataSize,'single');
    hdr.pixdim(2:4) = hdr.pixdim(2:4)./scalingFactor';
    hdr.sform44 = hdr.sform44*diag([1./scalingFactor 1]);
    hdr.qform44 = hdr.sform44*diag([1./scalingFactor 1]);


  
  
%read the warped fields
warpFields = mlrImageReadNifti(inversewarpfile);

scaledCoords = repmat([scalingFactor 1]',1,size(coords,2)).*(coords-.5) + .5;

%partition coordinates betwen those that are inside and outside the known warp fields
innerCoordsIndices = scaledCoords(1,:)>=1 & scaledCoords(1,:)<=dataSize(1) &...
                     scaledCoords(2,:)>=1 & scaledCoords(2,:)<=dataSize(2) &...
                     scaledCoords(3,:)>=1 & scaledCoords(3,:)<=dataSize(3);
                     
%interpolate fields
fieldVectors = zeros(4,size(coords,2));
innerCoords = scaledCoords(:,innerCoordsIndices);
%this might only work for isotropic voxels
fieldVectors(1,innerCoordsIndices) = interp3(warpFields(:,:,:,1),innerCoords(2,:),innerCoords(1,:),innerCoords(3,:),'*cubic'); %need to switch 1st and 2nd dimensions of the volume to call interp3
fieldVectors(2,innerCoordsIndices) = interp3(warpFields(:,:,:,2),innerCoords(2,:),innerCoords(1,:),innerCoords(3,:),'*cubic');
fieldVectors(3,innerCoordsIndices) = interp3(warpFields(:,:,:,3),innerCoords(2,:),innerCoords(1,:),innerCoords(3,:),'*cubic');


%the field values are presumably in millimeters, so needs to be divided by the voxel size t
fieldVectors = diag([coordsVoxelSize 1])\fieldVectors;
warpedCoords = coords-fieldVectors;

% %remove NaNs
% warpedCoords = warpedCoords(:,any(~isnan(warpedCoords),1));



