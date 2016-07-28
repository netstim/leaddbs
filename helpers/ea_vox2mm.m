function mm=ea_vox2mm(vox, transform)
% converts voxel-coordinates to mm-coordinates

if ischar(transform)
    transform = spm_get_space(transform);
end

vox = round(vox);
mm = [vox, ones(size(vox,1),1)] * transform';
mm(:,4) = [];
