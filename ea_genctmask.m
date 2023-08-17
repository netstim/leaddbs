function ea_genctmask(options)

masks_path = fileparts(options.subj.recon.rawCTMask);
ea_mkdir(masks_path);

if ~isfolder(masks_path); mkdir(masks_path); end

copyfile([ea_space,'brainmask.nii.gz'], fullfile(masks_path,'brainmask.nii.gz'));
gunzip(fullfile(masks_path,'brainmask.nii.gz'));
ea_delete(fullfile(masks_path,'brainmask.nii.gz'));

ea_apply_normalization_tofile(options, fullfile(masks_path,'brainmask.nii'), fullfile(masks_path,'wbrainmask.nii'), 1, 0);
ea_delete(fullfile(masks_path,'brainmask.nii'));

coreg_log = loadjson(options.subj.coreg.log.method);

if contains(coreg_log.method.CT, 'ANTs')
    transform_file_name = [options.subj.coreg.transform.CT.inverseBaseName 'ants.mat'];
elseif contains(coreg_log.method.CT, 'BRAINSFit')
    transform_file_name = [options.subj.coreg.transform.CT.inverseBaseName 'brainsfit.mat'];
elseif contains(coreg_log.method.CT, 'FLIRT')
    transform_file_name = [options.subj.coreg.transform.CT.inverseBaseName 'flirt.mat'];
end

ea_apply_coregistration(options.subj.preproc.anat.postop.CT, fullfile(masks_path,'wbrainmask.nii'), options.subj.recon.rawCTMask, ...
    transform_file_name,'nn'); % nn interpolation

movefile(fullfile(masks_path,'wbrainmask.nii'), options.subj.recon.anchorNativeMask);
