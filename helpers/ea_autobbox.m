function [bbox, BW] = ea_autobbox(image, margin)
% Calculate the minimum bounding box and binary mask of the nifti image

if nargin < 2
    margin = 5;	% add margin to the calculated bounding box, unit in voxel
end

nii = load_nii(image);
img = uint8(nii.img);

th = ea_otsuthresh(img);

BW = img > max(img(:)) * th;

[i,j,k] = ind2sub(size(BW), find(BW));

if strcmp(image(end-2:end), '.gz')
    gunzip(image);
    transform = spm_get_space(image(1:end-3));
    delete(image(1:end-3));
else
    transform = spm_get_space(image);
end

bbox = round([ea_vox2mm([min(i)-margin, min(j)-margin, min(k)-margin], transform);
              ea_vox2mm([max(i)+margin, max(j)+margin, max(k)+margin], transform)]);
