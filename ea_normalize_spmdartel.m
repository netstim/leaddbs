function varargout=ea_normalize_spmdartel(options)
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
% The procedure used here uses the SPM DARTEL approach to map a patient's
% brain to MNI space directly. Unlike the usual DARTEL-approach, which is
% usually used for group studies, here, DARTEL is used for a pairwise
% co-registration between patient anatomy and MNI template. It has been
% shown that DARTEL also performs superior to many other normalization approaches 
% also  in a pair-wise setting e.g. in 
%   Klein, A., et al. (2009). Evaluation of 14 nonlinear deformation algorithms 
%   applied to human brain MRI registration. NeuroImage, 46(3), 786?802. 
%   doi:10.1016/j.neuroimage.2008.12.037
%
% Since a high resolution is needed for accurate DBS localizations, this
% function applies DARTEL to an output resolution of 0.5 mm isotropic. This
% makes the procedure quite slow.

% The function uses some code snippets written by Ged Ridgway.
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn


if ischar(options) % return name of method.
    if strcmp(SPM('ver'),'SPM12')
        varargout{1}='SPM12 DARTEL nonlinear [MR/CT]';
    elseif strcmp(SPM('ver'),'SPM8')
        varargout{1}='SPM8 DARTEL nonlinear [MR/CT]';
    end
    varargout{2}={'SPM8','SPM12'};
    return
end
if ~exist([options.earoot,'templates',filesep,'TPM.nii'],'file')
   ea_generate_tpm; % will generate a hd_template out of SPM's TPM.nii
    
end

segmentresolution=0.5; % resolution of the DARTEL-Warps. Setting this value to larger values will generate the usual DARTEL-Workflow.



usecombined=0; % if set, eauto will try to fuse coronar and transversal images before normalizing them.
usesegmentnew=0;
costfuns={'nmi','mi','ecc','ncc'};

if exist([options.root,options.prefs.patientdir,filesep,options.prefs.tranii_unnormalized,'.gz'],'file')
    try
        gunzip([options.root,options.prefs.patientdir,filesep,options.prefs.tranii_unnormalized,'.gz']);
    end
    try
        gunzip([options.root,options.prefs.patientdir,filesep,options.prefs.cornii_unnormalized,'.gz']);
    end
    try
        gunzip([options.root,options.prefs.patientdir,filesep,options.prefs.sagnii_unnormalized,'.gz']);
    end
    
    
    try
        gunzip([options.root,options.prefs.patientdir,filesep,options.prefs.prenii_unnormalized,'.gz']);
    end
end

% check if backup files exist, if not backup

if ~exist([options.root,options.prefs.patientdir,filesep,'backup_',options.prefs.tranii_unnormalized],'file')
    try    copyfile([options.root,options.prefs.patientdir,filesep,options.prefs.tranii_unnormalized],[options.root,options.prefs.patientdir,filesep,'backup_',options.prefs.tranii_unnormalized]); end
else
    try copyfile([options.root,options.prefs.patientdir,filesep,'backup_',options.prefs.tranii_unnormalized],[options.root,options.prefs.patientdir,filesep,options.prefs.tranii_unnormalized]); end
end
if ~exist([options.root,options.prefs.patientdir,filesep,'backup_',options.prefs.cornii_unnormalized],'file')
    try    copyfile([options.root,options.prefs.patientdir,filesep,options.prefs.cornii_unnormalized],[options.root,options.prefs.patientdir,filesep,'backup_',options.prefs.cornii_unnormalized]); end
else
    try    copyfile([options.root,options.prefs.patientdir,filesep,'backup_',options.prefs.cornii_unnormalized],[options.root,options.prefs.patientdir,filesep,options.prefs.cornii_unnormalized]); end
end
if ~exist([options.root,options.prefs.patientdir,filesep,'backup_',options.prefs.sagnii_unnormalized],'file')
    try    copyfile([options.root,options.prefs.patientdir,filesep,options.prefs.sagnii_unnormalized],[options.root,options.prefs.patientdir,filesep,'backup_',options.prefs.sagnii_unnormalized]); end
else
    try    copyfile([options.root,options.prefs.patientdir,filesep,'backup_',options.prefs.sagnii_unnormalized],[options.root,options.prefs.patientdir,filesep,options.prefs.sagnii_unnormalized]); end
end


% First, do the coreg part:

    ea_coreg(options,options.prefs.normalize.coreg);




% now dartel-import the preoperative version.

%try



    disp('Segmenting preoperative version (Import to DARTEL-space)');
    
switch(spm('ver'))
    case 'SPM8'
        spmsegroot=[fileparts(which('spm')),filesep,'toolbox',filesep,'Seg',filesep];
    case 'SPM12'
        spmsegroot=[fileparts(which('spm')),filesep,'tpm',filesep];
end
    load([options.earoot,'ext_libs',filesep,'segment',filesep,'segjob']);
    job.channel.vols{1}=[options.root,options.prefs.patientdir,filesep,options.prefs.prenii_unnormalized];
    for tpm=1:6
        job.tissue(tpm).tpm=fullfile([spmsegroot,'TPM.nii,',num2str(tpm)]);
        if tpm<4
        job.tissue(tpm).native=[0,1];
        else
            job.tissue(tpm).native=[0,0];
        end
    end
    job.resolution=segmentresolution; 
    if strcmp(spm('ver'),'SPM12');
        job.warp.fwhm=0;
        job.warp.cleanup=1;
        job.warp.reg=[0 0.001 0.5 0.05 0.2];
        spm_preproc_run(job);
    else
        
            ea_spm_preproc_run(job); % exactly the same as the SPM version ("New Segment" in SPM8 / "Segment" in SPM12) but with an increase in resolution to 0.5 mm iso.
    end
    
    
    
    disp('*** Segmentation of preoperative MRI worked.');

    
    
    % check if darteltemplate is available, if not generate one
    
    if exist([options.earoot,filesep,'templates',filesep,'dartel',filesep,'dartelmni_6.nii'],'file')
        % There is a DARTEL-Template. Check if it will match:
        Vt=spm_vol([options.earoot,filesep,'templates',filesep,'dartel',filesep,'dartelmni_6.nii']);
        Vp=spm_vol([options.root,filesep,options.patientname,filesep,'rc1pre_tra.nii']);
        if ~isequal(Vp.dim,Vt(1).dim) || ~isequal(Vp.mat,Vt(1).mat) % Dartel template not matching. -> create matching one.
            ea_create_mni_darteltemplate([options.root,filesep,options.patientname,filesep,'rc1pre_tra.nii']);
        end
        
    else % no dartel template present. -> Create matching dartel templates from highres version.
        ea_create_mni_darteltemplate([options.root,filesep,options.patientname,filesep,'rc1pre_tra.nii']);
        
    end
    
    %
    
    


% Normalize to MNI using DARTEL.
matlabbatch{1}.spm.tools.dartel.warp1.images = {
                                                {[options.root,options.prefs.patientdir,filesep,'rc1',options.prefs.prenii_unnormalized,',1']}
                                                {[options.root,options.prefs.patientdir,filesep,'rc2',options.prefs.prenii_unnormalized,',1']}
                                                {[options.root,options.prefs.patientdir,filesep,'rc3',options.prefs.prenii_unnormalized,',1']}
                                                }';
matlabbatch{1}.spm.tools.dartel.warp1.settings.rform = 0;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(1).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(1).rparam = [4 2 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(1).K = 0;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(1).template = {[options.earoot,'templates',filesep,'dartel',filesep,'dartelmni_6.nii']};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).rparam = [2 1 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).K = 0;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).template = {[options.earoot,'templates',filesep,'dartel',filesep,'dartelmni_5.nii']};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).rparam = [1 0.5 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).K = 1;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).template = {[options.earoot,'templates',filesep,'dartel',filesep,'dartelmni_4.nii']};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).rparam = [0.5 0.25 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).K = 2;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).template = {[options.earoot,'templates',filesep,'dartel',filesep,'dartelmni_3.nii']};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).rparam = [0.25 0.125 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).K = 4;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).template = {[options.earoot,'templates',filesep,'dartel',filesep,'dartelmni_2.nii']};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).rparam = [0.25 0.125 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).K = 6;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).template = {[options.earoot,'templates',filesep,'dartel',filesep,'dartelmni_1.nii']};
matlabbatch{1}.spm.tools.dartel.warp1.settings.optim.lmreg = 0.01;
matlabbatch{1}.spm.tools.dartel.warp1.settings.optim.cyc = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.optim.its = 3;
jobs{1}=matlabbatch;
%try
    cfg_util('run',jobs);
    disp('*** Dartel coregistration of preoperative version worked.');
%catch
%    ea_error('*** Dartel coregistration failed.');
%end
clear matlabbatch jobs;

% export normalization parameters:
for inverse=0:1
    if inverse
        addstr='_inv';
    else
        addstr='';
    end
    matlabbatch{1}.spm.util.defs.comp{1}.dartel.flowfield = {[options.root,options.prefs.patientdir,filesep,'u_rc1',options.prefs.prenii_unnormalized]};
    matlabbatch{1}.spm.util.defs.comp{1}.dartel.times = [1-inverse 0+inverse];
    matlabbatch{1}.spm.util.defs.comp{1}.dartel.K = 6;
    matlabbatch{1}.spm.util.defs.ofname = ['ea',addstr,'_normparams'];
    matlabbatch{1}.spm.util.defs.fnames = '';
    matlabbatch{1}.spm.util.defs.savedir.saveusr = {[options.root,options.prefs.patientdir,filesep]};
    matlabbatch{1}.spm.util.defs.interp = 1;
    jobs{1}=matlabbatch;

    cfg_util('run',jobs);
    disp('*** Exported normalization parameters to y_ea_normparams.nii');
    clear matlabbatch jobs;
end


ea_apply_normalization(options)
