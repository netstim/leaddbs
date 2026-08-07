function srcIdx = ea_unified_perm_nnindexmap(fromSpace, toSpace)
% Build a nearest-neighbor voxel correspondence between two nifti grids,
% without touching disk / SPM reslice. Used to project a training object's
% saved permutation maps (thousands of value-rows sharing one fixed
% fromSpace/toSpace pair) onto a different test object's grid cheaply: the
% expensive part (which source voxel is nearest to which destination voxel)
% is computed once here, then applied to any number of value-vectors via a
% single indexing operation (permvals(valid) = row(srcIdx(valid))).
%
% Training and test grids are NOT guaranteed to match, even at the same
% resolution and even both nominally in MNI space: this codebase builds
% each analysis's grid as the union bounding box of that specific cohort's
% own VTA/efield extents (see ea_unifiedmapping_exportefieldmap.m), so two
% different patient cohorts generally produce different .dim/.mat even at
% identical resolution. At different resolutions there is no 1:1 index
% correspondence at all (e.g. one 1mm voxel roughly spans an 8-voxel block
% of a 0.5mm grid), so a geometric (world-mm) lookup is unavoidable there.
%
% Nearest-neighbor (not trilinear/smoothing) is used deliberately: Rperm
% rows are per-voxel correlation coefficients/p-values, not a continuous
% physical field, so blending values across a voxel boundary isn't
% meaningful -- every output voxel should equal an actual observed value
% from the source grid.
%
% fromSpace, toSpace - nii-like structs with .dim (1x3) and .mat (4x4),
%                       e.g. entries of obj.results.sweetspotmapping.space.
%
% srcIdx - numel(toSpace grid) x 1. srcIdx(j) is the linear voxel index
%          into fromSpace's flattened image nearest to destination voxel
%          j, or NaN if voxel j falls outside fromSpace's grid entirely.

dimTo = toSpace.dim(1:3);

% Fast path: grids are exactly identical (same dim AND same affine, not
% merely same resolution) -- the correct correspondence is the identity
% mapping, no geometric computation needed.
if isequal(fromSpace.dim(1:3), dimTo) && isequal(fromSpace.mat, toSpace.mat)
    srcIdx = (1:prod(dimTo))';
    return
end

[X,Y,Z] = ndgrid(1:dimTo(1), 1:dimTo(2), 1:dimTo(3));
voxTo = [X(:)'; Y(:)'; Z(:)'; ones(1,numel(X))]; % 4 x Nvoxels, 1-based voxel coords

worldXYZ = toSpace.mat * voxTo; % voxel -> world mm

voxFrom = fromSpace.mat \ worldXYZ; % world mm -> source voxel (continuous)
voxFrom = round(voxFrom(1:3,:)); % nearest neighbor

dimFrom = fromSpace.dim(1:3);
valid = voxFrom(1,:)>=1 & voxFrom(1,:)<=dimFrom(1) & ...
        voxFrom(2,:)>=1 & voxFrom(2,:)<=dimFrom(2) & ...
        voxFrom(3,:)>=1 & voxFrom(3,:)<=dimFrom(3);

srcIdx = nan(numel(X),1);
srcIdx(valid) = sub2ind(dimFrom, voxFrom(1,valid), voxFrom(2,valid), voxFrom(3,valid));
