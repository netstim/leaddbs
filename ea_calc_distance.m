function varargout=ea_calc_distance(varargin)
% Example:
% [normfactor,dist_mm,dist_vox]=calc_distance(priordist,trajectory,mat,nii)
% priordist is the distance in real-life in mm
% trajectory is the vector on which the electrodes lie
% mat is the transformation matrix that has been used to change the
% distances.
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

priordist=varargin{1};
trajvector=varargin{2};
mat=varargin{3};

newtrajvector=trajvector*mat; % i.e. the vector between 0,0,0 and some point in direction of trajvector.
% stretch factor in trajectory as found by mat.
normfactor=norm(newtrajvector)/norm(trajvector); % scalar value denoting the relationship between distances defined by trajvector in native and mni space.
% distance in millimeters
dist_mm=priordist*normfactor;

varargout{2}=dist_mm;
varargout{1}=normfactor;

% calculate voxel distances if nargin > 3.
if nargin>3
    nii=varargin{4};
    XYZ_vx = ea_mm2vox([0,0,0;trajvector], nii);

    trajvox=diff(XYZ_vx)';

    worldtovoxfactor=norm(trajvox)/norm(trajvector);
    dist_vox=priordist*worldtovoxfactor*normfactor;
    varargout{3}=dist_vox;
end
