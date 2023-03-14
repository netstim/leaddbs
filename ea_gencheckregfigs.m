function ea_gencheckregfigs(options, type)
fprintf('\nExporting coregistration check images to %scheckreg...\n', [options.subj.coregDir,filesep]);

if ~exist('type', 'var')
    type = {'coreg', 'norm'};
elseif ischar(type)
    type = {type};
end

if ismember('coreg', type)
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

    % Get paths of coregistered image
    coregImage = [struct2cell(preop); struct2cell(postop)];

    % Get paths of output figures
    if isempty(struct2cell(preop)) && isempty(struct2cell(postop))
        ea_cprintf('CmdWinWarnings', 'No coregistration performed, checkcoreg shipped!\n');
        return;
    elseif ~isempty(struct2cell(preop)) && ~isempty(struct2cell(postop))
        checkregFigure = [struct2cell(options.subj.coreg.checkreg.preop); struct2cell(options.subj.coreg.checkreg.postop)];
    elseif isempty(struct2cell(preop)) % Only one pre-op modality available
        checkregFigure = struct2cell(options.subj.coreg.checkreg.postop);
    elseif isempty(struct2cell(postop)) % Only one pre-op modality available
        checkregFigure = struct2cell(options.subj.coreg.checkreg.preop);
    end

    % Generate checkreg figures
    for i=1:length(coregImage)
        if isfile(coregImage{i})
            ea_gencheckregpair(coregImage{i}, anchorImage, checkregFigure{i});
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
        ea_gencheckregpair(options.subj.brainshift.anat.moving, anchorImage, options.subj.brainshift.checkreg.standard);
    end

    % Generate checkreg figure for brain shift corrected image
    if isfile(options.subj.brainshift.anat.scrf)
        ea_gencheckregpair(options.subj.brainshift.anat.scrf, anchorImage, options.subj.brainshift.checkreg.scrf);
    end
end

if ismember('norm', type)
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

fprintf('Done.\n');
