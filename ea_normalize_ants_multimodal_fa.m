function varargout=ea_normalize_ants_multimodal_fa(options)
% This is a function that normalizes both a copy of transversal and coronar
% images into MNI-space. The goal was to make the procedure both robust and
% automatic, but still, it must be said that normalization results should
% be taken with much care because all reconstruction results heavily depend
% on these results. Normalization of DBS-MR-images is especially
% problematic since usually, the field of view doesn't cover the whole
% brain (to reduce SAR-levels during acquisition) and since electrode
% artifacts can impair the normalization process. Therefore, normalization
% might be best archieved with other tools that have specialized on
% normalization of such image data.
%
% The procedure used here uses the ANTs Syn approach to map a patient's
% brain to MNI space directly. 
% __________________________________________________________________________________
% Copyright (C) 2015 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn


if ischar(options) % return name of method.
    varargout{1}='Advanced Normalization Tools (ANTs) SyN - multimodal, include FA';
    varargout{2}={'SPM8','SPM12'};
    return
end


directory=[options.root,options.patientname,filesep];



% check for presence of FA map
if ~exist([directory,options.prefs.fa2anat],'file')
    if ~exist([directory,options.prefs.fa],'file')
        if ~exist([directory,options.prefs.dti],'file')
           ea_error('Please provide a DTI image (e.g. dti.nii) in subject folder.'); 
        end
        ea_isolate_fa(directory,options);
    end
    switch options.coregmr.method
        case 1 % SPM
            ea_docoreg_spm([directory,options.prefs.fa],[directory,options.prefs.prenii_unnormalized],'nmi',1)
            movefile([directory,'r',options.prefs.fa],[directory,options.prefs.fa2anat]);
        case 2 % ANTs
            ea_ants([directory,options.prefs.prenii_unnormalized],...
                [directory,options.prefs.fa],...
                [directory,options.prefs.fa2anat],0);
        case 3 % BRAINSFit
            ea_brainsfit([directory,options.prefs.prenii_unnormalized],...
                [directory,options.prefs.fa],...
                [directory,options.prefs.fa2anat],0);
        case 4 % Hybrid SPM -> ANTs
            ea_docoreg_spm([directory,options.prefs.fa],[directory,options.prefs.prenii_unnormalized],'nmi',0)
            ea_ants([directory,options.prefs.prenii_unnormalized],...
                [directory,options.prefs.fa],...
                [directory,options.prefs.fa2anat],0);
        case 5 % Hybrid SPM -> Brainsfit
            ea_docoreg_spm([directory,options.prefs.fa],[directory,options.prefs.prenii_unnormalized],'nmi',0)
            ea_brainsfit([directory,options.prefs.prenii_unnormalized],...
                [directory,options.prefs.fa],...
                [directory,options.prefs.fa2anat],0);
    end
end

% First, do the coreg part:

%ea_coregmr(options,options.prefs.normalize.coreg);

keyboard

ea_ants_nonlinear({[options.earoot,'templates',filesep,'mni_hires.nii'],...
    [options.earoot,'templates',filesep,'mni_hires_t1.nii'],...
    [options.earoot,'templates',filesep,'mni_hires_pd.nii']
    [options.earoot,'templates',filesep,'mni_hires_fa.nii']},...
    {[directory,options.prefs.prenii_unnormalized],...
    [directory,options.prefs.prenii_unnormalized],...
    [directory,options.prefs.prenii_unnormalized]
    [directory,options.prefs.fa2anat]},...
    [directory,options.prefs.prenii]);


ea_apply_normalization(options)
