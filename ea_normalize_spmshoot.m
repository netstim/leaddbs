function varargout=ea_normalize_spmshoot(options)
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
    varargout{1}='SPM12 SHOOT (Ashburner 2011)';
    varargout{2}=1; % dummy output
    varargout{3}=0; % hassettings.
    varargout{4}=1; % is multispectral
    return
end

disp('Segmenting preoperative version (Import to DARTEL-space)');
preopImages = struct2cell(options.subj.coreg.anat.preop);
ea_newseg(preopImages, 1, 1);
disp('Segmentation of preoperative MRI done.');

[directory, preopAnchorName] = fileparts(preopImages{1});
directory = [directory, filesep];

% Check if SHOOT template is available
if exist([ea_space(options,'dartel'),'shootmni_6.nii'],'file')
    % There is a DARTEL-Template. Check if it will match:
    Vt=spm_vol([ea_space(options,'dartel'),'shootmni_6.nii']);
    Vp=spm_vol([directory, 'rc1', preopAnchorName, '.nii']);
    if ~isequal(Vp.dim,Vt(1).dim) || ~isequal(Vp.mat,Vt(1).mat) % Dartel template not matching. -> create matching one.
        ea_create_tpm_darteltemplate;
    end
else % no dartel template present. -> Create matching dartel templates from highres version.
    keyboard
    ea_create_tpm_darteltemplate;
end

% Forward
matlabbatch{1}.spm.tools.shoot.warp1.images = {
    {[directory,'rc1',preopAnchorName,'.nii,1']}
    {[directory,'rc2',preopAnchorName,'.nii,1']}
    {[directory,'rc3',preopAnchorName,'.nii,1']}}';
matlabbatch{1}.spm.tools.shoot.warp1.templates = {[ea_space(options,'dartel'),'shootmni_1.nii']
    [ea_space(options,'dartel'),'shootmni_2.nii']
    [ea_space(options,'dartel'),'shootmni_3.nii']
    [ea_space(options,'dartel'),'shootmni_4.nii']
    [ea_space(options,'dartel'),'shootmni_5.nii']
    [ea_space(options,'dartel'),'shootmni_6.nii']};

spm_jobman('run',{matlabbatch});
disp('*** Shoot coregistration of preoperative version worked.');
clear matlabbatch;

movefile([directory,'y_rc1', preopAnchorName, '_Template.nii'], [directory,'y_ea_normparams.nii']);

% Inverse
matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def = {[directory,'y_ea_normparams.nii']};
matlabbatch{1}.spm.util.defs.comp{1}.inv.space = {[directory, preopAnchorName, '.nii']};
matlabbatch{1}.spm.util.defs.out{1}.savedef.ofname = 'ea_inv_normparams.nii';
matlabbatch{1}.spm.util.defs.out{1}.savedef.savedir.saveusr = {directory};
spm_jobman('run',{matlabbatch});
disp('*** Exported normalization parameters to y_ea_inv_normparams.nii');
clear matlabbatch;

% Delete rc* files and u_rc1* file
ea_delete([directory, 'rc*', preopAnchorName, '.nii']);
ea_delete([directory, '*_rc1', preopAnchorName, '_Template.nii']);

% Rename Segmentations (c1, c2, c3)
mod = replace(options.subj.AnchorModality, textBoundary('start') + alphanumericsPattern + "_", "");
movefile([directory, 'c1', preopAnchorName, '.nii'], setBIDSEntity(preopImages{1}, 'mod', mod, 'label', 'GM', 'suffix', 'mask'));
movefile([directory, 'c2', preopAnchorName, '.nii'], setBIDSEntity(preopImages{1}, 'mod', mod, 'label', 'WM', 'suffix', 'mask'));
movefile([directory, 'c3', preopAnchorName, '.nii'], setBIDSEntity(preopImages{1}, 'mod', mod, 'label', 'CSF', 'suffix', 'mask'));

% Rename deformation fields
ea_mkdir(fileparts(options.subj.norm.transform.forwardBaseName));
movefile([directory, 'y_ea_normparams.nii'], [options.subj.norm.transform.forwardBaseName, 'spm.nii']);
movefile([directory, 'y_ea_inv_normparams.nii'], [options.subj.norm.transform.inverseBaseName, 'spm.nii']);

ea_apply_normalization(options)

% add methods dump:
[scit, lcit] = ea_getspacedefcit;
cits={
    'Ashburner, J., & Friston, K. J. (2005). Unified segmentation., 26(3), 839?851. http://doi.org/10.1016/j.neuroimage.2005.02.018'
    'Ashburner, J., & Friston, K. J. (2011). Diffeomorphic registration using geodesic shooting and Gauss?Newton optimisation. NeuroImage, 55(3), 954?967. http://doi.org/10.1016/j.neuroimage.2010.12.049'
    'Horn, A., & Kuehn, A. A. (2015). Lead-DBS: a toolbox for deep brain stimulation electrode localizations and visualizations. NeuroImage, 107, 127?135. http://doi.org/10.1016/j.neuroimage.2014.12.002'};

if ~isempty(lcit)
    cits = [cits; {lcit}];
end

modality = regexp(preopImages, '(?<=_)[^\W_]+(?=\.nii(\.gz)?$)', 'match', 'once');

ea_methods(options,['Pre- (and post-) operative acquisitions were spatially normalized into ',ea_getspace,' space ',scit,' based on preoperative acquisition(s) (',strjoin(modality, ', '),') using a'...
    ' diffeomorphic registration algorithm using geodesic shooting and Gauss-Neuwton optimisation (SHOOT) as implemented in SPM12 (Ashburner 2011; www.fil.ion.ucl.ac.uk/spm/software/).',...
    ' SHOOT registration was performed by directly registering tissue segmentations of preoperative acquisitions (obtained using the unified Segmentation approach as implemented in SPM12 (Ashburner 2005)',...
    ' to a SHOOT template created from tissue priors defined by the MNI (ICBM 152 Nonlinear asymmetric 2009b atlas; http://nist.mni.mcgill.ca/?p=904)',...
    ' supplied within Lead-DBS software (Horn 2015; www.lead-dbs.org).'],cits);
