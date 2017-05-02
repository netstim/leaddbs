function [bbox, BW] = ea_autobbox(image, margin)
% Calculate the minimum bounding box and binary mask of the nifti image

if nargin < 2
    margin = 3;	% add margin to the calculated bounding box, unit in voxel
end

nii = load_nii(image);
img = uint8(nii.img);

th = otsuthresh(img);

BW = img > max(img(:)) * th;

[i,j,k] = ind2sub(size(BW), find(BW));

bbox = round([ea_vox2mm([min(i)-margin, min(j)-margin, min(k)-margin], image);
              ea_vox2mm([max(i)+margin, max(j)+margin, max(k)+margin], image)]);
