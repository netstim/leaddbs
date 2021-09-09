function varargout=ea_normalize_spmdartel(options)
% This is a function that normalizes both a copy of transversal and coronal
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
    varargout{1}='SPM12 DARTEL (Ashburner 2007)';
    varargout{2}=1; % dummy output
    varargout{3}=0; % hassettings.
    varargout{4}=1; % is multispectral
    return
end


directory = [options.root,options.patientname,filesep];

if isfield(options.prefs, 'tranii_unnormalized')
    if exist([directory,options.prefs.tranii_unnormalized,'.gz'],'file')
        try
            gunzip([directory,options.prefs.tranii_unnormalized,'.gz']);
        end
        try
            gunzip([directory,options.prefs.cornii_unnormalized,'.gz']);
        end
        try
            gunzip([directory,options.prefs.sagnii_unnormalized,'.gz']);
        end
        try
            gunzip([directory,options.prefs.prenii_unnormalized,'.gz']);
        end
    end
end



% now dartel-import the preoperative version.
disp('Segmenting preoperative version (Import to DARTEL-space)');
ea_newseg_pt(options,1,1);
disp('Segmentation of preoperative MRI done.');

% check if darteltemplate is available, if not generate one
if exist([ea_space(options,'dartel'),'dartelmni_6.nii'],'file')
    % There is a DARTEL-Template. Check if it will match:
    Vt=spm_vol([ea_space(options,'dartel'),'dartelmni_6.nii']);
    Vp=spm_vol([directory,'rc1',options.prefs.prenii_unnormalized]);
    if ~isequal(Vp.dim,Vt(1).dim) || ~isequal(Vp.mat,Vt(1).mat) % Dartel template not matching. -> create matching one.
        ea_create_tpm_darteltemplate; %([directory,'rc1',options.prefs.prenii_unnormalized]);
    end
else % no dartel template present. -> Create matching dartel templates from highres version.
    ea_create_tpm_darteltemplate; %([directory,'rc1',options.prefs.prenii_unnormalized]);
end

% Normalize to MNI using DARTEL.
matlabbatch{1}.spm.tools.dartel.warp1.images = {
                                                {[directory,'rc1',options.prefs.prenii_unnormalized,',1']}
                                                {[directory,'rc2',options.prefs.prenii_unnormalized,',1']}
                                                {[directory,'rc3',options.prefs.prenii_unnormalized,',1']}
                                                }';
matlabbatch{1}.spm.tools.dartel.warp1.settings.rform = 0;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(1).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(1).rparam = [4 2 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(1).K = 0;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(1).template = {[ea_space(options,'dartel'),'dartelmni_1.nii']};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).rparam = [2 1 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).K = 0;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(2).template = {[ea_space(options,'dartel'),'dartelmni_2.nii']};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).rparam = [1 0.5 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).K = 1;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(3).template = {[ea_space(options,'dartel'),'dartelmni_3.nii']};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).rparam = [0.5 0.25 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).K = 2;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(4).template = {[ea_space(options,'dartel'),'dartelmni_4.nii']};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).rparam = [0.25 0.125 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).K = 4;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(5).template = {[ea_space(options,'dartel'),'dartelmni_5.nii']};
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).its = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).rparam = [0.25 0.125 1e-06];
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).K = 6;
matlabbatch{1}.spm.tools.dartel.warp1.settings.param(6).template = {[ea_space(options,'dartel'),'dartelmni_6.nii']};
matlabbatch{1}.spm.tools.dartel.warp1.settings.optim.lmreg = 0.01;
matlabbatch{1}.spm.tools.dartel.warp1.settings.optim.cyc = 3;
matlabbatch{1}.spm.tools.dartel.warp1.settings.optim.its = 3;
jobs{1}=matlabbatch;

spm_jobman('run',jobs);
disp('*** Dartel coregistration of preoperative version worked.');

clear matlabbatch jobs;

% Export normalization parameters:
% backward
switch spm('ver')
    case 'SPM8'
        ea_error('SPM8 is not supported anymore in this version of Lead-DBS');
    case 'SPM12'
        matlabbatch{1}.spm.util.defs.comp{1}.dartel.flowfield = {[directory,'u_rc1',options.prefs.prenii_unnormalized]};
        matlabbatch{1}.spm.util.defs.comp{1}.dartel.times = [1 0];
        matlabbatch{1}.spm.util.defs.comp{1}.dartel.K = 6;
        matlabbatch{1}.spm.util.defs.comp{1}.dartel.template = {''};
        matlabbatch{1}.spm.util.defs.out{1}.savedef.ofname = ['ea_normparams'];
        matlabbatch{1}.spm.util.defs.out{1}.savedef.savedir.saveusr = {directory};
end
jobs{1}=matlabbatch;

spm_jobman('run',jobs);
disp('*** Exported normalization parameters to y_ea_normparams.nii');
clear matlabbatch jobs;

% forward (inverse)
switch spm('ver')
    case 'SPM8'
        ea_error('SPM8 is not supported anymore in this version of Lead-DBS');
    case 'SPM12'
        matlabbatch{1}.spm.util.defs.comp{1}.dartel.flowfield = {[directory,'u_rc1',options.prefs.prenii_unnormalized]};
        matlabbatch{1}.spm.util.defs.comp{1}.dartel.times = [0 1];
        matlabbatch{1}.spm.util.defs.comp{1}.dartel.K = 6;
        matlabbatch{1}.spm.util.defs.comp{1}.dartel.template = {''};
        matlabbatch{1}.spm.util.defs.comp{2}.id.space = {[directory,options.prefs.prenii_unnormalized]};
        matlabbatch{1}.spm.util.defs.out{1}.savedef.ofname = 'ea_inv_normparams';
        matlabbatch{1}.spm.util.defs.out{1}.savedef.savedir.saveusr = {directory};
end
jobs{1}=matlabbatch;

spm_jobman('run',jobs);
disp('*** Exported normalization parameters to y_ea_inv_normparams.nii');
clear matlabbatch jobs;

% delete([directory,'u_rc1',options.prefs.prenii_unnormalized]);

ea_apply_normalization(options)

%% add methods dump:
[scit,lcit]=ea_getspacedefcit;
cits={
    'Ashburner, J., & Friston, K. J. (2005). Unified segmentation., 26(3), 839?851. http://doi.org/10.1016/j.neuroimage.2005.02.018'
    'Ashburner, J. (2007). A fast diffeomorphic image registration algorithm, 38(1), 95?113. http://doi.org/10.1016/j.neuroimage.2007.07.007'
    'Horn, A., & Kuehn, A. A. (2015). Lead-DBS: a toolbox for deep brain stimulation electrode localizations and visualizations. NeuroImage, 107, 127?135. http://doi.org/10.1016/j.neuroimage.2014.12.002'
    };
if ~isempty(lcit)
    cits=[cits;{lcit}];
end
[~,anatpresent]=ea_assignpretra(options);

ea_methods(options,['Pre- (and post-) operative acquisitions were spatially normalized into ',ea_getspace,' space ',scit,' based on preoperative acquisition(s) (',ea_cell2strlist(anatpresent),') using a'...
    ' fast diffeomorphic image registration algorithm (DARTEL) as implemented in ',spm('ver'),' (Ashburner 2007; www.fil.ion.ucl.ac.uk/spm/software/).',...
    ' DARTEL registration was performed by directly registering tissue segmentations of preoperative acquisitions (obtained using the unified Segmentation approach as implemented in ',spm('ver'),' (Ashburner 2005)',...
    ' to a DARTEL template created from tissue priors defined by the MNI (ICBM 152 Nonlinear asymmetric 2009b atlas; http://nist.mni.mcgill.ca/?p=904)',...
    ' supplied within Lead-DBS software (Horn 2015; www.lead-dbs.org).',...
    ],...
    cits);
