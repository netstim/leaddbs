function mm=vox2mm(vox, transform)

% converts voxel-coordinates to mm-coordinates

if ischar(transform)
    v2m = spm_get_space(transform);
else
    v2m = transform;
end

vox = round(vox);
mm = [vox, ones(size(vox,1),1)] * v2m';
mm(:,4) = [];
