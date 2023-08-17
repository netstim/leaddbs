function varargout = ea_normalize_schoenecker(options)
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
% The procedure used here uses the ANTs Syn approach to map a patient's
% brain to MNI space directly.
% __________________________________________________________________________________
% Copyright (C) 2015 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn


if ischar(options) % return name of method.
    varargout{1}='Three-step affine normalization (ANTs; Schonecker 2009)';
    varargout{2}=1; % dummy output
    varargout{3}=1; % hassettings.
    varargout{4}=1; % is multispectral
    return
end

usefa = options.prefs.machine.normsettings.ants_usefa;
usebrainmask=0;

cnt=1;

% TODO: Take care of FA
spacedef = options.bids.spacedef;
if usefa && spacedef.hasfa % first put in FA since least important (if both an FA template and an fa2anat file is available)
    if exist(options.prefs.fa2anat,'file') % recheck if now is present.
        disp('Including FA information for white-matter normalization.');
        template{cnt} = [ea_space(options),'fa.nii'];
        moving{cnt} = [bprfx,options.prefs.fa2anat];
        weights(cnt) = 0.5;
        metrics{cnt} = 'MI';
        cnt = cnt+1;
    end
end

if options.prefs.machine.normsettings.schoenecker_movim==1 % Based on pre-op images
    disp(['Pre-op ', strjoin(fieldnames(options.subj.coreg.anat.preop), ', '), ' images included for normalization']);
    imagePresent = flip(struct2cell(options.subj.coreg.anat.preop)); % Flip the order so anchor will be the last one
elseif options.prefs.machine.normsettings.schoenecker_movim==2 % Based on post-op images
    switch options.subj.postopModality
        case 'MRI'
            disp('Post-op MRI included for normalization');
            imagePresent = options.subj.coreg.anat.postop.ax_MRI;
        case 'CT'
            disp('Post-op CT included for normalization');
            imagePresent = options.subj.coreg.anat.postop.CT;
    end
end

% The convergence criterion for the multivariate scenario is a slave to the
% last metric you pass on the ANTs command line.
for i=1:length(imagePresent)
    % Set template image
    template{cnt} = ea_matchTemplate(imagePresent{i}, spacedef);

    % Segment and mask image when needed.
    if usebrainmask
        moving{cnt} = ea_maskimg(imagePresent{i});
    else
        moving{cnt} = imagePresent{i};
    end

    % Set weights and metrics
    weights(cnt) = 1.25;
    metrics{cnt} = 'MI';
    cnt = cnt+1;
end

ea_ants_schoenecker(template, moving, options.subj.norm.anat.preop.(options.subj.AnchorModality), weights, metrics);

% Move transformation file
[directory, fileName] = fileparts(options.subj.norm.anat.preop.(options.subj.AnchorModality));
ea_mkdir(fileparts(options.subj.norm.transform.forwardBaseName));
movefile([directory, filesep, fileName, '0GenericAffine.mat'], [options.subj.norm.transform.forwardBaseName, 'ants.mat']);
movefile([directory, filesep, fileName, 'Inverse0GenericAffine.mat'], [options.subj.norm.transform.inverseBaseName, 'ants.mat']);

ea_apply_normalization(options);

% Add methods dump
[scit, lcit] = ea_getspacedefcit;
cits={
    'Avants, B. B., Epstein, C. L., Grossman, M., & Gee, J. C. (2008). Symmetric diffeomorphic image registration with cross-correlation: evaluating automated labeling of elderly and neurodegenerative brain. Medical Image Analysis, 12(1), 26?41. http://doi.org/10.1016/j.media.2007.06.004'
    'Schonecker, T., Kupsch, A., Kuhn, A. A., Schneider, G. H., & Hoffmann, K. T. (2009). Automated optimization of subcortical cerebral MR imaging-atlas coregistration for improved postoperative electrode localization in deep brain stimulation. AJNR American Journal of Neuroradiology, 30(10), 1914?1921. http://doi.org/10.3174/ajnr.A1741'};

if ~isempty(lcit)
    cits=[cits;{lcit}];
end

switch options.prefs.machine.normsettings.schoenecker_movim
    case 1 % pre-op image used
        methodVerbose = ' Unlike in the original publication, the three-step registration was estimated based on the preoperative acquisition and applied to the (co-registered) postoperative acquisition.';
    case 2 % post-op image used
        methodVerbose = '';
end

modality = regexp(imagePresent, '(?<=_)[^\W_]+(?=\.nii(\.gz)?$)', 'match', 'once');
ea_methods(options,['Pre- (and post-) operative acquisitions were spatially normalized into ',ea_getspace,' space ',scit,' based on preoperative acquisition(s) (',strjoin(modality, ', '),') using a',...
    ' three-step linear affine registration as implemented in Advanced Normalization Tools (Avants 2008; http://stnava.github.io/ANTs/). This linear registration protocol follows the approach described in Schonecker et al. 2009.',...
    ' After a brief rigid-body transform, linear deformation into template space was achieved in three stages: 1. Whole brain affine registration 2. Affine registration that solved the cost-function inside a subcortical mask only and',...
    ' 3. Affine registration that focused on a region defined by a small basal ganglia mask. In the original publication, this approach was validated for use in DBS and yielded a',...
    ' RMS error of 1.29 mm.', methodVerbose], cits);
