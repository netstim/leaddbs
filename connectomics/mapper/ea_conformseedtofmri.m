function seed = ea_conformseedtofmri(dataset,seed)
td = tempdir;
tmpref = [td,'tmpspace.nii'];
tmpseed = [td,'tmpseed.nii'];

% Write out ref into temp file
dataset.vol.space.fname = tmpref;
ea_write_nii(dataset.vol.space);

% Write out seed into temp file
seed.fname = tmpseed;
ea_write_nii(seed);

% Conform space by coregistration using SPM with trilinear interpolation
options.coregmr.method='SPM';
ea_coreg2images(options,tmpseed,tmpref,tmpseed,[],[],[],1);

% Load resliced seed file and clean up
seed=ea_load_nii(tmpseed);
delete(tmpref);
delete(tmpseed);
