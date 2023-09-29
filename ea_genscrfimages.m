function ea_genscrfimages(subj, type)
% Generated brainshift corrected postop images under coregistration/anat
% and normalization/anat.
arguments
    subj {mustBeA(subj, {'struct', 'char', 'string'})}
    type {mustBeMember(type, {'coreg', 'norm', 'all'})} = 'all'
end

if isstruct(subj)
    % Input is BIDS struct
    options.subj = subj;
    spacedef = ea_getspacedef;
    options.primarytemplate = fullfile(ea_space, [spacedef.templates{1}, '.nii']);
elseif isfolder(subj)
    % Input is patient derivatives folder
    options = ea_getptopts(subj);
    subj = options.subj;
end

postopFields = fieldnames(subj.postopAnat);
for i=1:length(postopFields)
    % Generate brainshift corrected coreg image
    if ismember(type, {'coreg', 'all'})
        ea_cprintf('*Comments', '\nGenerating brainshift corrected coreg image: %s ...\n', postopFields{i});
        ea_ants_apply_transforms(struct, ...
            subj.postopAnat.(postopFields{i}).coreg, ... % Moving
            subj.postopAnat.(postopFields{i}).coregScrf, ... % Output
            0, ... % Forward transform
            subj.preopAnat.(subj.AnchorModality).coreg, ... % Reference
            subj.brainshift.transform.instore, ... % Transform
            'Linear');
    
        if strcmp(subj.postopModality, 'CT')
            % Generate brainshift corrected coreg image (tone-mapped CT)
            ea_cprintf('*Comments', '\nGenerating brainshift corrected coreg image: tone-mapped %s ...\n', postopFields{i});
            ea_ants_apply_transforms(struct, ...
                subj.postopAnat.(postopFields{i}).coregTonemap, ... % Moving
                subj.postopAnat.(postopFields{i}).coregTonemapScrf, ... % Output
                0, ... % Forward transform
                subj.preopAnat.(subj.AnchorModality).coreg, ... % Reference
                subj.brainshift.transform.instore, ... % Transform
                'Linear');
        end
    end

    % Generate brainshift corrected norm image
    if ismember(type, {'norm', 'all'})
        ea_cprintf('*Comments', '\nGenerating brainshift corrected norm image: %s ...\n', postopFields{i});
        ea_apply_normalization_tofile(options, ...
            subj.postopAnat.(postopFields{i}).coregScrf, ... % Moving
            subj.postopAnat.(postopFields{i}).normScrf, ... % Output
            0, ... % Forward transform
            1, ... % Linear interpolation
            options.primarytemplate); % Reference
    
        if strcmp(subj.postopModality, 'CT')
            % Generate brainshift corrected norm image (tone-mapped CT)
            ea_cprintf('*Comments', '\nGenerating brainshift corrected norm image: tone-mapped %s ...\n', postopFields{i});
            ea_apply_normalization_tofile(options, ...
                subj.postopAnat.(postopFields{i}).coregTonemapScrf, ... % Moving
                subj.postopAnat.(postopFields{i}).normTonemapScrf, ... % Output
                0, ... % Forward transform
                1, ... % Linear interpolation
                options.primarytemplate); % Reference
        end
    end
end
