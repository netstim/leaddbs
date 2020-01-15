function vox=ea_mm2vox(mm, transform)
% converts mm-coordinates to voxel-coordinates
% coords need to be row vector: N*3

if ischar(transform)
    transform = spm_get_space(transform);
end

vox = [mm, ones(size(mm,1),1)] / transform';
vox(:,4) = [];
% vox = round(vox);
