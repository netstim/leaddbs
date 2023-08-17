function seed = ea_conformseedtofmri(dataset,seed,interp)
% Use SPM coreg to conform space for VTA.
    
% Can solve the all nan/zero problem for some very small VTA. The course is
% that rigid body reslicing provided by spm_reslice is not enough. 12 DOF
% (affine) or 9 DOF (traditional) transformation might be needed.

if ~isstruct(dataset) && isfile(dataset)
    dataset = load(dataset);
end

if ~isstruct(seed) && isfile(seed)
    seed = ea_load_nii(seed);
end

% Set default interp method to trilinear
% Can choose from trilinear, nearestneighbour, sinc or spline
if ~exist('interp', 'var')
    interp = 'trilinear';
end
guid=ea_generate_uuid;
td = tempdir;
tmpref = [td,'tmpspace_',guid,'.nii'];
tmpseed = [td,'tmpseed_',guid,'.nii'];

% Write out ref into temp file
dataset.vol.space.fname = tmpref;
dataset.vol.space.dt(1) = 16;
ea_write_nii(dataset.vol.space);

% Write out seed into temp file
seed.fname = tmpseed;
seed.dt(1) = 16;
ea_write_nii(seed);

% Conform space using FSL or SPM reslicing with trilinear interpolation
prefs=ea_prefs;
switch prefs.lcm.vat2fmrimethod
    case 'fsl' % default since v2.5.3
        ea_fsl_reslice(tmpseed, tmpref, tmpseed, interp);
    case 'spm' % older method - can lead to issues with very small vtas.
        ea_conformspaceto(tmpref, tmpseed, 1); % use interp = 1 for trilinear.
end

% Load resliced seed file and clean up
seed=ea_load_nii(tmpseed);
delete(tmpref);
delete(tmpseed);
