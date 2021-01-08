function seed = ea_conformseedtofmri(dataset,seed,interp)
% Use SPM coreg to conform space for VTA.
    
% Can solve the all nan/zero problem for some very small VTA. The course is
% that rigid body reslicing provided by spm_reslice is not enough. 12 DOF
% (affine) or 9 DOF (traditional) transformation might be needed.

if ~isstruct(dataset) && isfile(dataset)
    load(dataset, 'dataset');
end

if ~isstruct(seed) && isfile(seed)
    seed = ea_load_nii(seed);
end

% Use trilinear interpolation by default.
% Set to 0 to use nearest neighbour
if ~exist('interp', 'var')
    interp = 1;
end

td = tempdir;
tmpref = [td,'tmpspace.nii'];
tmpseed = [td,'tmpseed.nii'];

% Write out ref into temp file
dataset.vol.space.fname = tmpref;
dataset.vol.space.dt = [16,0];
ea_write_nii(dataset.vol.space);

% Write out seed into temp file
seed.fname = tmpseed;
seed.dt = [16,0];
ea_write_nii(seed);

% Conform space by coregistration using SPM with trilinear interpolation
options.coregmr.method='SPM';
ea_coreg2images(options,tmpseed,tmpref,tmpseed,[],[],[],interp);

% Load resliced seed file and clean up
seed=ea_load_nii(tmpseed);
delete(tmpref);
delete(tmpseed);
