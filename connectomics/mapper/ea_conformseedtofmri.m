function s = ea_conformseedtofmri(dataset,s)
td = tempdir;
ref = [td,'tmpspace.nii'];
seed = [td,'tmpseed.nii'];

% Write out ref into temp file
dataset.vol.space.fname = ref;
ea_write_nii(dataset.vol.space);

% Write out seed into temp file
s.fname = seed;
ea_write_nii(s);

% Conform space by coregistration using SPM with trilinear interpolation
options.coregmr.method='SPM';
ea_coreg2images(options,seed,ref,seed,[],[],[],1);

% Load resliced seed file and clean up
s=ea_load_nii(seed);
delete(ref);
delete(seed);
