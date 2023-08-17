function varargout = ea_normalize_ants(options)
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
    varargout{1} = 'ANTs (Avants 2008)';
    varargout{2} = 1; % dummy output
    varargout{3} = 1; % hassettings.
    varargout{4} = 1; % is multispectral
    return
end

usefa = options.prefs.machine.normsettings.ants_usefa;
usebrainmask=0;

cnt=1;

% TODO: Take care of FA
spacedef = options.bids.spacedef;
if usefa && spacedef.hasfa % first put in FA since least important (if both an FA template and an fa2anat file is available)
    if exist([options.prefs.fa2anat],'file') % recheck if now is present.
        disp('Including FA information for white-matter normalization.');
        template{cnt} = [ea_space(options),'fa.nii'];
        moving{cnt} = [bprfx,options.prefs.fa2anat];
        weights(cnt) = 0.5;
        cnt = cnt+1;
    elseif exist([options.prefs.fa],'file') % recheck if now is present.
        disp('Including FA information for white-matter normalization (weight = 0.5).');
        ea_coregimages(options,[bprfx,options.prefs.fa],[anatpresent{1}],[bprfx,options.prefs.fa2anat],{},0,[],1);
        template{cnt} = [ea_space(options),'fa.nii'];
        moving{cnt} = [bprfx,options.prefs.fa2anat];
        weights(cnt) = 0.5;
        cnt = cnt+1;
    end
end

disp(['Pre-op ', strjoin(fieldnames(options.subj.coreg.anat.preop), ', '), ' images included for normalization']);
imagePresent = flip(struct2cell(options.subj.coreg.anat.preop)); % Flip the order so anchor will be the last one

% The convergence criterion for the multivariate scenario is a slave to the
% last metric you pass on the ANTs command line.
for i = 1:length(imagePresent)
    % Set template image
    template{cnt} = ea_matchTemplate(imagePresent{i}, spacedef);

    % Segment and mask image when needed.
    if usebrainmask
        moving{cnt} = ea_maskimg(imagePresent{i});
    else
        moving{cnt} = imagePresent{i};
    end

    % Set weights
    weights(cnt) = 1.25;
    cnt = cnt+1;
end

% Add segmentations in case present (under coregistration/segmentations)
segmentations = ea_regexpdir([options.subj.subjDir, filesep, 'segmentations'], '.*\.nii(\.gz)?$', 0);
if ~isempty(segmentations)
    for i=1:length(segmentations)
        [~, segName] = ea_niifileparts(segmentations{i});
        segTemplate = ea_regexpdir([ea_space, 'segmentations'], [segName, '\.nii(\.gz)?$'], 0);
        if ~isempty(segTemplate) % Found matching template in template space segmentations folder
            disp(['Including ', segName, ' segmentation for segment based assistance (weight = 20).']);
            moving = [segmentations{i}, moving]; % append to front (since last one is convergence critical)
            template = [segTemplate{1}, template];
            weights = [20, weights]; % set weight to 20 - DO NOT CHANGE THIS VALUE BELOW 3. IF VALUE IS CHANGED, SEGMENTATIONS WILL BE CONSIDERED SLABS IN ea_ants_nonlinear ~line 63 - would need to be changed there, as well.
        end
    end
end

output_image = options.subj.norm.anat.preop.(options.subj.AnchorModality);
output_transform_prefix = fullfile(fileparts(options.subj.norm.transform.forwardBaseName), 'antsout');

ea_mkdir(fileparts(output_image));
ea_mkdir(fileparts(output_transform_prefix));

ea_ants_nonlinear(template, moving, output_image, weights, output_transform_prefix, options);

% Move transformation file
movefile([output_transform_prefix 'Composite.nii.gz'], [options.subj.norm.transform.forwardBaseName, 'ants.nii.gz']);
movefile([output_transform_prefix 'InverseComposite.nii.gz'], [options.subj.norm.transform.inverseBaseName, 'ants.nii.gz']);

ea_apply_normalization(options);

% Add methods dump
modality = regexp(imagePresent, '(?<=_)[^\W_]+(?=\.nii(\.gz)?$)', 'match', 'once');
if options.prefs.machine.normsettings.ants_scrf
    [scit, lcit] = ea_getspacedefcit;
    cits={
        'Avants, B. B., Epstein, C. L., Grossman, M., & Gee, J. C. (2008). Symmetric diffeomorphic image registration with cross-correlation: evaluating automated labeling of elderly and neurodegenerative brain. Medical Image Analysis, 12(1), 26?41. http://doi.org/10.1016/j.media.2007.06.004'
        'Schoenecker, T., Kupsch, A., Kuehn, A. A., Schneider, G.-H., & Hoffmann, K. T. (2009). Automated Optimization of Subcortical Cerebral MR Imaging-Atlas Coregistration for Improved Postoperative Electrode Localization in Deep Brain Stimulation. AJNR Am J Neuroradiol, 30(10), 1914?1921. http://doi.org/10.3174/ajnr.A1741'};
    if ~isempty(lcit)
        cits = [cits;{lcit}];
    end

    ea_methods(options,['Pre- (and post-) operative acquisitions were spatially normalized into ',ea_getspace,' space ',scit,' based on preoperative acquisition(s) (',strjoin(modality, ', '),') using the'...
        ' SyN registration approach as implemented in Advanced Normalization Tools (Avants 2008; http://stnava.github.io/ANTs/).',...
        ' Nonlinear deformation into template space was achieved in five stages: After two linear (rigid followed by affine) steps, ',...
        ' A nonlinear (whole brain) SyN-registration stage was followed by two nonlinear SyN-registrations that consecutively focused on the area of interest ',...
        ' as defined by subcortical masks in Schoenecker 2008.'],cits);
else
    [scit, lcit] = ea_getspacedefcit;
    cits = {'Avants, B. B., Epstein, C. L., Grossman, M., & Gee, J. C. (2008). Symmetric diffeomorphic image registration with cross-correlation: evaluating automated labeling of elderly and neurodegenerative brain. Medical Image Analysis, 12(1), 26?41. http://doi.org/10.1016/j.media.2007.06.004'};
    if ~isempty(lcit)
        cits = [cits;{lcit}];
    end

    ea_methods(options,['Pre- (and post-) operative acquisitions were spatially normalized into ',ea_getspace,' space ',scit,' based on preoperative acquisition(s) (',strjoin(modality, ', '),') using the'...
        ' SyN registration approach as implemented in Advanced Normalization Tools (Avants 2008; http://stnava.github.io/ANTs/).',...
        ' Nonlinear deformation into template space was achieved in three stages: After two linear (rigid followed by affine) steps, ',...
        ' a nonlinear (whole brain) SyN registration stage was added.'],cits);
end

