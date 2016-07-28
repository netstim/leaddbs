function vox=ea_mm2vox(mm, transform)
% converts mm-coordinates to voxel-coordinates

if ischar(transform)
    transform = spm_get_space(transform);
end
    
transform = inv(transform);
vox = [mm, ones(size(mm,1),1)] * transform';
vox(:,4) = [];
vox = round(vox);
