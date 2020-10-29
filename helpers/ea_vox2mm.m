function mm=ea_vox2mm(vox, transform, base)
% converts voxel-coordinates to mm-coordinates
% coords need to be row vector: N*3

if ischar(transform)
    if ~exist('base', 'var')
        base = 1; % Use one-based indexing by default
    end
    transform = ea_get_affine(transform, base);
end

% vox = round(vox);
mm = [vox, ones(size(vox,1),1)] * transform';
mm(:,4) = [];
