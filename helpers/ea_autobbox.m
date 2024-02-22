function [bbox, BW] = ea_autobbox(input, margin)
% Calculate the minimum bounding box and binary mask of the nifti image

if nargin < 2
    margin = 5;	% add margin to the calculated bounding box, unit in voxel
end

nii = ea_load_nii(input);
img = nii.img;

th = multithresh(img);

BW = img > th;

[i,j,k] = ind2sub(size(BW), find(BW));

transform = ea_get_affine(input);

bbox = round([ea_vox2mm([min(i)-margin, min(j)-margin, min(k)-margin], transform);
              ea_vox2mm([max(i)+margin, max(j)+margin, max(k)+margin], transform)]);

bbox = sort(bbox);
