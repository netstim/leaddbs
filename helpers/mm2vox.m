function vox=mm2vox(mm, transform)

% converts mm-coordinates to voxel-coordinates

if ischar(transform)
    v2m = spm_get_space(transform);
else
    v2m = transform;
end
    
m2v = inv(v2m);
vox = [mm, ones(size(mm,1),1)] * m2v';
vox(:,4) = [];
vox = round(vox);
