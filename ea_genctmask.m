function ea_genctmask(options)
directory=[options.root, options.patientname, filesep];
copyfile([ea_space,'brainmask.nii.gz'],[directory,'brainmask.nii.gz']);
gunzip([directory,'brainmask.nii.gz']);
ea_delete([directory,'brainmask.nii.gz']);

                             % options,           from,                 to,                          directory,useinverse,interp,refim)
ea_apply_normalization_tofile(options, {[directory,'brainmask.nii']}, {[directory,'wbrainmask.nii']}, directory, 1, 0);
ea_delete([directory,'brainmask.nii']);

load([directory,'ea_coregctmethod_applied.mat']) % determine last used coregmethod
switch coregct_method_applied{end}
    case 'ea_coregctmri_brainsfit'
        suffix='_brainsfit.h5';
    case {'ea_coregctmri_ants', 'ea_coregctmri_ants_refine'}
        coregs=dir([directory,ea_stripext(options.prefs.prenii_unnormalized),'2',ea_stripext(options.prefs.rawctnii_unnormalized),'_ants*.mat']);
        suffix=strrep(coregs(end).name,[ea_stripext(options.prefs.prenii_unnormalized),'2',ea_stripext(options.prefs.rawctnii_unnormalized)],'');
    case 'ea_coregctmri_fsl'
        coregs=dir([directory,ea_stripext(options.prefs.prenii_unnormalized),'2',ea_stripext(options.prefs.rawctnii_unnormalized),'_flirt*.mat']);
        suffix=strrep(coregs(end).name,[ea_stripext(options.prefs.prenii_unnormalized),'2',ea_stripext(options.prefs.rawctnii_unnormalized)],'');
end
ea_apply_coregistration([directory,options.prefs.rawctnii_unnormalized], [directory,'wbrainmask.nii'], [directory,'ct_mask.nii'], ...
    [directory,ea_stripext(options.prefs.prenii_unnormalized),'2',ea_stripext(options.prefs.rawctnii_unnormalized),suffix],'nn'); % nn interpolation
movefile([directory,'wbrainmask.nii'],[directory,'rct_mask.nii']);
