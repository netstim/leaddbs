function [bbox, BW] = ea_autobbox(image, margin)
% Calculate the minimum bounding box and binary mask of the nifti image

if nargin < 2
    margin = 5;	% add margin to the calculated bounding box, unit in voxel
end

img = uint8(spm_read_vols(spm_vol(image)));

th = ea_otsuthresh(img);

BW = img > max(img(:)) * th;

[i,j,k] = ind2sub(size(BW), find(BW));

transform = ea_get_affine(image);

bbox = round([ea_vox2mm([min(i)-margin, min(j)-margin, min(k)-margin], transform);
              ea_vox2mm([max(i)+margin, max(j)+margin, max(k)+margin], transform)]);
