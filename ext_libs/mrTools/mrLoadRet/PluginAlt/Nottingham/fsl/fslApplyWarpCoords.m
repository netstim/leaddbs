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

function warpedCoords = fslApplyWarpCoords(coords,coordsVoxelSize,warpResolution, warpCoefFilename, tempFilename, hdr, verbose)


if ieNotDefined('verbose')
   verbose = true;
end

 scalingFactor = round(hdr.pixdim(2:4)'./warpResolution);
    dataSize = hdr.dim(2:4)'.*scalingFactor;
    data = ones(dataSize,'single');
    hdr.pixdim(2:4) = hdr.pixdim(2:4)./scalingFactor';
    hdr.sform44 = hdr.sform44*diag([1./scalingFactor 1]);
    hdr.qform44 = hdr.sform44*diag([1./scalingFactor 1]);

if ~exist(tempFilename,'file') % only needs to be done once..
   
    
    cbiWriteNifti(tempFilename, data, hdr,[],[],[],verbose);
    clear data
    
    if ispc
        FNIRTUTILS=ea_path_helper([ea_getearoot,'ext_libs',filesep,'fsl',filesep,'fnirtfileutils','.exe']);
    else
        FNIRTUTILS=ea_path_helper([ea_getearoot,'ext_libs',filesep,'fsl',filesep,'fnirtfileutils','.', computer('arch')]);
    end
    
    command =  sprintf(' --in=%s --ref=%s --out=%s  --withaff', warpCoefFilename, tempFilename, tempFilename);
    if verbose
        fprintf('(fslApplyWarpCoords) Computing FNIRT warp fields at a resolution of %s mm:\n',mat2str(hdr.pixdim(2:4)));
        disp(['  ' command])
    end;
    lcmd=[FNIRTUTILS,command];
    if ~ispc
        system(['bash -c "', lcmd, '"']);
    else
        system(lcmd);
    end
end
  
  
%read the warped fields
warpFields = mlrImageReadNifti(tempFilename);

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


if any(~innerCoordsIndices)   %extrapolate fields outside FNIRT volume (ad-hoc method)
  
  outerCoordsIndices = find(~innerCoordsIndices);
  outerCoords = scaledCoords(:,outerCoordsIndices);
  
  %first estimate the linear component of the registration
  %using ordinary least squares
  %this will fill the warp fields outside the outer box
  [volumeXCoords,volumeYCoords,volumeZCoords] = ndgrid(1:dataSize(1),1:dataSize(2),1:dataSize(3));
  volumeCoords = [volumeXCoords(:) volumeYCoords(:) volumeZCoords(:) ones(prod(dataSize),1)]; %create regression matrix and add intercept column
  clear volumeXCoords volumeYCoords volumeZCoords
  %precompute the normal equations
  normalEquations = (volumeCoords'*volumeCoords)\volumeCoords';
  clear volumeCoords
  %compute the parameters to go from coordinates to linear fields
  betaDx = normalEquations*reshape(warpFields(:,:,:,1),prod(dataSize),1);
  betaDy = normalEquations*reshape(warpFields(:,:,:,2),prod(dataSize),1);
  betaDz = normalEquations*reshape(warpFields(:,:,:,3),prod(dataSize),1);
  clear normalEquations

  %we also need  a buffer zone between the non-linear fields and the linear fields around
  
  % This zone will be defined by an inner box exactly enclosing the volume
  % and an outer box at some separation distance around the FNIRT input volume
  separation =ceil(2*squeeze(max(max(max(max(warpFields)))))); %the buffer zone is twice the size of the largest field value

  %compute linear field values outside the input FNIRT volume
  fieldVectors(1,outerCoordsIndices) = outerCoords'*betaDx;
  fieldVectors(2,outerCoordsIndices) = outerCoords'*betaDy;
  fieldVectors(3,outerCoordsIndices) = outerCoords'*betaDz;
  
  %identify points that are between the outer and inner boxes
  betweenCoordsIndices = find(outerCoords(1,:)>=1-separation & outerCoords(1,:)<=dataSize(1)+separation &...
                     outerCoords(2,:)>=1-separation & outerCoords(2,:)<=dataSize(2)+separation &...
                     outerCoords(3,:)>=1-separation & outerCoords(3,:)<=dataSize(3)+separation);
  
  betweenCoords = outerCoords(:,betweenCoordsIndices);                 
      
  %project betweenCoords onto the inner (parallel for points outside one side 
  %and onto the edges/corners for points outside 2/3 sides 
  innerBoxCoords = betweenCoords;
  innerBoxCoords(1,:) = min(max(innerBoxCoords(1,:),1),dataSize(1));
  innerBoxCoords(2,:) = min(max(innerBoxCoords(2,:),1),dataSize(2));
  innerBoxCoords(3,:) = min(max(innerBoxCoords(3,:),1),dataSize(3));
  
  %project point onto outer box along the segment defined by the coordinate and its inner projection
  % fix separation distance between inner and outer box, which produces round edges
  
  projectionVectors = betweenCoords(1:3,:) - innerBoxCoords(1:3,:);
  innerDistance = sqrt(sum(projectionVectors.^2,1));
  %get rid of points that might be between the round and sharp edges of the outer box
  betweenCoords(:,innerDistance>separation)=[];
  betweenCoordsIndices(innerDistance>separation)=[];
  innerBoxCoords(:,innerDistance>separation)=[];
  projectionVectors(:,innerDistance>separation)=[];
  innerDistance(innerDistance>separation)=[];
  %project 
  outerBoxCoords = ones(size(innerBoxCoords));
  outerBoxCoords(1:3,:) = innerBoxCoords(1:3,:)+projectionVectors./repmat(innerDistance,3,1)*separation;
    
  %get interpolated values for projected points on inner box 
  %the values at this box are interpolated from the actual warp fields at the limit of the known volume
  innerBoxDxValues = interp3(warpFields(:,:,:,1),innerBoxCoords(2,:),innerBoxCoords(1,:),innerBoxCoords(3,:),'*cubic'); %need to switch 1st and 2nd dimensions of the volume to call interp3
  innerBoxDyValues = interp3(warpFields(:,:,:,2),innerBoxCoords(2,:),innerBoxCoords(1,:),innerBoxCoords(3,:),'*cubic');
  innerBoxDzValues = interp3(warpFields(:,:,:,3),innerBoxCoords(2,:),innerBoxCoords(1,:),innerBoxCoords(3,:),'*cubic');
  
  %compute values for projected points on outer box (from estimated linear transform)
  outerBoxDxValues = betaDx'*outerBoxCoords;
  outerBoxDyValues = betaDy'*outerBoxCoords;
  outerBoxDzValues = betaDz'*outerBoxCoords;
  
  %now intrepolate values between the inner and outer boxes (between the inner and outer projections)
  fieldVectors(1,outerCoordsIndices(betweenCoordsIndices)) = innerBoxDxValues+ (outerBoxDxValues-innerBoxDxValues).*innerDistance/separation;
  fieldVectors(2,outerCoordsIndices(betweenCoordsIndices)) = innerBoxDyValues+ (outerBoxDyValues-innerBoxDyValues).*innerDistance/separation;
  fieldVectors(3,outerCoordsIndices(betweenCoordsIndices)) = innerBoxDzValues+ (outerBoxDzValues-innerBoxDzValues).*innerDistance/separation;

end

%the field values are presumably in millimeters, so needs to be divided by the voxel size t
fieldVectors = diag([coordsVoxelSize 1])\fieldVectors;
warpedCoords = coords-fieldVectors;

% %remove NaNs
% warpedCoords = warpedCoords(:,any(~isnan(warpedCoords),1));



