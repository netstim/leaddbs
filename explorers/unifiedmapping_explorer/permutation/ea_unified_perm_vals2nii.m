function fnames = ea_unified_perm_vals2nii(space, vals, outdir, basename)
% Write a voxelwise stat map (one column vector per side) into nifti
% file(s), using the grid geometry in space{side} (e.g.
% obj.results.sweetspotmapping.space, or the .space field saved in a
% permtest results .mat file).
%
% space    - cell array of nii structs (need .dim/.mat/.img; .dt and .img
%            get overwritten here, so the source space{side}.img itself is
%            untouched by this call).
% vals     - cell array, one Nvoxels x 1 column vector per side, matching
%            numel(space{side}.img).
% outdir   - output folder (created if it doesn't exist).
% basename - filename prefix, e.g. 'sweetspot_Remp'. Side suffix ('_r'/'_l')
%            is appended automatically when there are 2 sides.
%
% Returns fnames, the cell array of written file paths.

ea_mkdir(outdir);

if numel(space) == 2
    sidesuffices = {'_r', '_l'};
else
    sidesuffices = {''};
end

fnames = cell(1, numel(space));
for side = 1:numel(space)
    nii = space{side};
    nii.dt(1) = 16; % float32, these are continuous stat values (R/p), not efield magnitudes
    nii.img(:) = vals{side};
    nii.fname = fullfile(outdir, [basename, sidesuffices{side}, '.nii']);
    ea_write_nii(nii);
    fnames{side} = nii.fname;
end
