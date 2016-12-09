function varargout=ea_normalize_ants_macaque(options)
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
    varargout{1}='Advanced Normalization Tools (ANTs) SyN';
    varargout{2}={'SPM8','SPM12'};
    return
end



% do a linear coregistration into mni space
if options.modality==1 %MR
    expdo=2:4;
elseif options.modality==2 % CT
    expdo=6;
end

cnt=1;
for export=expdo
    switch export
        case 2
            if exist([options.root,options.prefs.patientdir,filesep,options.prefs.tranii_unnormalized],'file')
                finas{cnt}=[options.root,options.prefs.patientdir,filesep,options.prefs.tranii_unnormalized];
                cnt=cnt+1;
            end
        case 3
            if exist([options.root,options.prefs.patientdir,filesep,options.prefs.cornii_unnormalized],'file')
                finas{cnt}=[options.root,options.prefs.patientdir,filesep,options.prefs.cornii_unnormalized];
                cnt=cnt+1;
            end
        case 4
            if exist([options.root,options.prefs.patientdir,filesep,options.prefs.sagnii_unnormalized],'file')
                finas{cnt}=[options.root,options.prefs.patientdir,filesep,options.prefs.sagnii_unnormalized];
                cnt=cnt+1;
            end
        case 6 % CT
            if exist([options.root,options.prefs.patientdir,filesep,options.prefs.ctnii_coregistered],'file')
                finas{cnt}=[options.root,options.prefs.patientdir,filesep,options.prefs.ctnii_coregistered];
                cnt=cnt+1;
            end
    end
end

matlabbatch{1}.spm.spatial.coreg.estimate.ref = {fullfile(options.earoot,'toolbox','macaque','templates','mni_hires_t2.nii,1')};
matlabbatch{1}.spm.spatial.coreg.estimate.source = {[options.root,options.prefs.patientdir,filesep,options.prefs.prenii_unnormalized,',1']};
matlabbatch{1}.spm.spatial.coreg.estimate.other = finas;
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
spm_jobman('run',{matlabbatch});

directory=[options.root,options.patientname,filesep];
ea_ants_nonlinear([options.earoot,'toolbox',filesep,'macaque',filesep,'templates',filesep,'mni_hires','.nii'],[directory,options.prefs.prenii_unnormalized],[directory,options.prefs.gprenii]);

ea_apply_normalization(options)
