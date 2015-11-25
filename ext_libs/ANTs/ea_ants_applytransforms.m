function ea_ants_applytransforms(options)
% Wrapper for antsApplyTransforms in terms of reapplying normalizations to
% pre- and postop imaging.

ea_libs_helper;

basedir = [fileparts(mfilename('fullpath')), filesep];

if ispc
    applyTransforms = [basedir, 'antsApplyTransforms.exe'];
elseif isunix
    applyTransforms = [basedir, 'antsApplyTransforms.', computer];
end

subdir=[options.root,options.patientname,filesep];

switch options.modality
    case 1 % MR
        fis{1}=[subdir,options.prefs.prenii_unnormalized];
        fis{2}=[subdir,options.prefs.tranii_unnormalized];
        fis{3}=[subdir,options.prefs.cornii_unnormalized];
        fis{4}=[subdir,options.prefs.sagnii_unnormalized];
        ofis{1}=[subdir,options.prefs.prenii];
        ofis{2}=[subdir,options.prefs.tranii];
        ofis{3}=[subdir,options.prefs.cornii];
        ofis{4}=[subdir,options.prefs.sagnii];
    case 2 % CT
        fis{1}=[subdir,options.prefs.prenii_unnormalized];
        fis{2}=[subdir,options.prefs.ctnii_coregistered];
        ofis{1}=[subdir,options.prefs.prenii];
        ofis{2}=[subdir,options.prefs.ctnii];
end

for fi=1:length(fis)
    [~,lprebase]=fileparts(options.prefs.prenii);
    system([applyTransforms,' --verbose 1' ...
        ' --dimensionality 3 --float 1' ...
        ' -i ',fis{fi} ...
        ' -o ',ofis{fi} ...
        ' -r ',[options.earoot,'templates',filesep,'mni_hires.nii']...
        ' -t ',[subdir,lprebase],'1Warp.nii.gz'...
        ' -t ',[subdir,lprebase],'0GenericAffine.mat']);
    
end