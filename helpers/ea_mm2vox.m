function vox=ea_mm2vox(mm, transform, base)
% converts mm-coordinates to voxel-coordinates
% coords need to be row vector: N*3

if ischar(transform)
    if ~exist('base', 'var')
        base = 1; % Use one-based indexing by default
    end
    transform = ea_get_affine(transform, base);
end

vox = [mm, ones(size(mm,1),1)] / transform';
vox(:,4) = [];
% vox = round(vox);
