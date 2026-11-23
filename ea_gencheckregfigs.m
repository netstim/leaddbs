function ea_gencheckregfigs(options, type)

if ~exist('type', 'var')
    type = {'coreg', 'norm'};
elseif ischar(type)
    type = {type};
end

if ismember('coreg', type)
    fprintf('\nExporting coregistration check images to %scheckreg...\n', [options.subj.coregDir,filesep]);
    % Make dir
    ea_mkdir(fullfile(options.subj.coregDir, 'checkreg'));

    % Get anchor Image path
    anchorImage = options.subj.coreg.anat.preop.(options.subj.AnchorModality);

    % Remove anchor image from subj.coreg.anat.preop struct
    preop = rmfield(options.subj.coreg.anat.preop, options.subj.AnchorModality);

    % Remove CT from subj.coreg.anat.postop struct, will use tone-mapped image
    if strcmp(options.subj.postopModality, 'CT')
        postop = rmfield(options.subj.coreg.anat.postop, 'CT');
    elseif strcmp(options.subj.postopModality, 'MRI')
        postop = options.subj.coreg.anat.postop;
    else
        postop = struct;
    end

    if isempty(fieldnames(preop)) && isempty(fieldnames(postop))
        ea_cprintf('CmdWinWarnings', 'No coregistration performed, checkcoreg shipped!\n');
        return;
    end

    % Pair each coregistered image with its checkreg figure by modality
    % name (not by position, which breaks as soon as one struct has an
    % entry, e.g. b0 or fa, that the other one lacks).
    sessions = {'preop', preop; 'postop', postop};
    for s = 1:size(sessions, 1)
        session = sessions{s, 1};
        images = sessions{s, 2};
        if ~isfield(options.subj.coreg.checkreg, session)
            continue;
        end
        figures = options.subj.coreg.checkreg.(session);
        modalities = fieldnames(images);
        for m = 1:numel(modalities)
            coregImage = images.(modalities{m});
            if ~isfield(figures, modalities{m})
                ea_cprintf('CmdWinWarnings', 'No checkreg figure path defined for %s (%s), skipped.\n', modalities{m}, session);
                continue;
            end
            if isfile(coregImage)
                ea_gencheckregpair(coregImage, anchorImage, figures.(modalities{m}));
            end
        end
    end
end

if ismember('norm', type)
    fprintf('\nExporting normalization check images to %scheckreg...\n', [options.subj.normDir,filesep]);
    % Make dir
    ea_mkdir(fullfile(options.subj.normDir, 'checkreg'));

    % Get template Image path
    templateImage = [ea_space, options.primarytemplate, '.nii'];

    % Remove CT from subj.norm.anat.postop struct, will use tone-mapped image
    if strcmp(options.subj.postopModality, 'CT')
        postop = rmfield(options.subj.norm.anat.postop, 'CT');
    elseif strcmp(options.subj.postopModality, 'MRI')
        postop = options.subj.norm.anat.postop;
    else
        postop = struct;
    end

    % Get paths of coregistered image
    normImage = [struct2cell(options.subj.norm.anat.preop); struct2cell(postop)];

    % Get paths of output figures
    if ~isempty(struct2cell(postop))
        checkregFigure = [struct2cell(options.subj.norm.checkreg.preop); struct2cell(options.subj.norm.checkreg.postop)];
    else
        checkregFigure = struct2cell(options.subj.norm.checkreg.preop);
    end

    % Generate checkreg figures
    for i=1:length(normImage)
        if isfile(normImage{i})
            ea_gencheckregpair(normImage{i}, templateImage, checkregFigure{i});
        end
    end
end

if ismember('brainshift', type)
    % Make dir
    ea_mkdir(fullfile(options.subj.brainshiftDir, 'checkreg'));

    % Get anchor Image path
    anchorImage = options.subj.brainshift.anat.anchor;

    % Generate checkreg figure for standard image
    if isfile(options.subj.brainshift.anat.moving)
        fprintf('\nExporting brainshift correction check image (standard) to %scheckreg...\n', [options.subj.brainshiftDir,filesep]);
        ea_gencheckregpair(options.subj.brainshift.anat.moving, anchorImage, options.subj.brainshift.checkreg.standard, preproc = 0);
    end

    % Generate checkreg figure for brain shift corrected image
    if isfile(options.subj.brainshift.anat.scrf)
        fprintf('\nExporting brainshift correction check image (scrf) to %scheckreg...\n', [options.subj.brainshiftDir,filesep]);
        ea_gencheckregpair(options.subj.brainshift.anat.scrf, anchorImage, options.subj.brainshift.checkreg.scrf, preproc = 0);
    end
end

fprintf('Done.\n');
