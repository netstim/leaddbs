function roi = ea_spherical_roi(fname,center,radius,crop,ref,bg)
% Create sphere ROI based on specified center and radius (both in mm)

% Write out NIfTI or not
if isempty(fname)
    writeoutNii = 0;
else
    writeoutNii = 1;
end

% Expand radius in case multiple centers specified
if size(center,1)>1
    if length(radius)==1
        radius = repmat(radius, 1, size(center,1));
    elseif size(center,1)~=length(radius)
        error('Length of centers doesn''t match length of radius!');
    end
end

% Crop the generate ROI image or not
if ~exist('crop','var')
    crop=1;
end

% Reference template image, use MNI t1 by default
if exist('ref','var')
    ref = ea_load_nii(ref);
else
    ref = ea_load_nii([ea_space,'t1.nii']);
end

% Preset background
if ~exist('bg','var')
    ref.img(:) = 0;
else
    if isscalar(bg)
        ref.img(:) = bg;
    elseif isequal(size(ref.img), size(bg))
        ref.img = bg;
    else
        error('Background should be either a scalar value or of the same size as the reference image!');
    end
end

voxsize = ref.voxsize;
dim = ref.dim;

for i=1:size(center,1)
    % mm to voxel conversion
    c = ea_mm2vox(center(i,:), ref.mat);
    r = radius(i)./voxsize;

    % Construct voxel grid for the sphere of cencter c and radius r
    bboxlim = [max([1 1 1; ceil(c-r)]); min([dim; floor(c+r)])];
    [xgrid, ygrid, zgrid] = meshgrid(bboxlim(1,1):bboxlim(2,1),...
                                     bboxlim(1,2):bboxlim(2,2),...
                                     bboxlim(1,3):bboxlim(2,3));

    % Flatten voxel grid to x, y and z subscripts
    xyz = [xgrid(:), ygrid(:), zgrid(:)];

    % Find voxels within the sphere
    xyz = xyz(sqrt(sum(((xyz - c) .* voxsize) .^ 2, 2)) <= radius(i), :);
    ref.img(sub2ind(dim, xyz(:,1), xyz(:,2), xyz(:,3))) = 1;
end

% Adapt ROI NIfTI structure
ref.dt(1) = 16;
ref.fname = fname;

roi = ref;

% Write out NIfTI
if writeoutNii
    ea_write_nii(ref);
    % Crop ROI image
    if crop
        ea_autocrop(fname);
    end
end
