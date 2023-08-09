function ea_createscrfmask(options)
% Create refining basal ganglia masks for brain shift correction

masks = {options.subj.brainshift.anat.secondstepmask
    options.subj.brainshift.anat.thirdstepmask};

if ~all(isfile(masks)) || options.overwriteapproved
    disp('Registering subcortical mask (Schoenecker 2008) to subject pre-op anatomy...');

    copyfile([ea_space(options,'subcortical'), 'secondstepmask.nii'], masks{1});
    copyfile([ea_space(options,'subcortical'), 'thirdstepmask.nii'], masks{2});

    anchor = options.subj.preopAnat.(options.subj.AnchorModality).coreg;

    for i=1:2
        % Warp masks
        if isfile(options.subj.preopAnat.(options.subj.AnchorModality).norm) % If normalization has been done already, use inverse warp
            ea_apply_normalization_tofile(options,masks(i),masks(i),1,0,anchor);
        else % If not, use a simple linear coregistration
            ea_coregimages(options, [ea_space, options.primarytemplate, '.nii'],...
                anchor, [fileparts(masks{i}), filesep, 'tmp.nii'], masks(i), 0);
            ea_delete([fileparts(masks{i}), filesep, 'tmp.nii']);
        end

        % Confrom space the brain shift anchor image
        ea_conformspaceto(options.subj.brainshift.anat.anchor, masks{i}, 0);
    end
end
