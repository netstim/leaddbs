function roi = ea_spherical_roi(fname,center,radius,crop,ref,bg)

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
if exist('template','var')
    ref = ea_load_nii(ref);
else
    ref = ea_load_nii([ea_space,'t1.nii']);
end

% Preset background
if ~exist('bg','var')
    ref.img(:) = 0;
else
    ref.img = bg;
end

voxsize = ref.voxsize;
for i=1:size(center,1)
    % mm to voxel conversion
    c = ea_mm2vox(center(i,:), ref.mat);
    r = radius(i);

    % Span along axes
    xspan = round(r/voxsize(1))*2 + 1;
    yspan = round(r/voxsize(2))*2 + 1;
    zspan = round(r/voxsize(3))*2 + 1;

    % Create grid, Contruct sphere within the grid
    [xgrid, ygrid, zgrid] = meshgrid(1:xspan,1:yspan,1:zspan);
    S = sqrt((xgrid-r/voxsize(1)).^2+(ygrid-r/voxsize(2)).^2+(zgrid-r/voxsize(3)).^2)<=r/mean(voxsize);

    % Relocate grid in the image space
    xgrid = xgrid + round(c(1)-r/voxsize(1)-1);
    ygrid = ygrid + round(c(2)-r/voxsize(2)-1);
    zgrid = zgrid + round(c(3)-r/voxsize(3)-1);

    % Fix grid outside of the image space
    xgrid(xgrid<1) = 1;
    xgrid(xgrid>size(ref.img,1)) = size(ref.img,1);
    ygrid(ygrid<1) = 1;
    ygrid(ygrid>size(ref.img,2)) = size(ref.img,2);
    zgrid(zgrid<1) = 1;
    zgrid(zgrid>size(ref.img,3)) = size(ref.img,3);

    % Convert grid to indices in the image space
    gridInd = sub2ind(size(ref.img), xgrid, ygrid, zgrid);

    % Find image indices within the sphere
    sphereInd = unique(gridInd(S));

    % Set sphere ROI in image
    ref.img(sphereInd) = 1;
end

% Adapt ROI NIfTI structure
ref.dt = [16,0];
ref.fname = fname;

roi = ref;

% Write out NIfTI
if writeoutNii
    ea_write_nii(ref);
    % Crop ROI image
    if crop
        ea_autocrop(fname)
    end
end
