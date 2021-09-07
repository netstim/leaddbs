function ea_gencheckregfigs(options, type)
fprintf('\nExporting coregistration check images to %scheckreg...\n', [options.root,options.patientname,filesep]);

if ~exist('type', 'var')
    type = {'coreg', 'norm'};
elseif ischar(type)
    type = {type};
end

if ismember('coreg', type)
    % Make dir
    ea_mkdir(fullfile(options.subj.subjDir, 'coregistration', 'checkreg'));

    % Get anchor Image path
    anchorImage = options.subj.coreg.anat.preop.(options.subj.AnchorModality);

    % Remove anchor image from subj.coreg.anat.preop struct
    preop = rmfield(options.subj.coreg.anat.preop, options.subj.AnchorModality);

    % Remove CT from subj.coreg.anat.postop struct
    if options.modality
        postop = rmfield(options.subj.coreg.anat.postop, 'CT');
    end

    % Get paths of coregistered image
    coregImage = [struct2cell(preop); struct2cell(postop)];

    % Get paths of output figures
    checkregFigure = [struct2cell(options.subj.coreg.checkreg.preop); struct2cell(options.subj.coreg.checkreg.postop)];

    % Generate checkreg figures
    for i=1:length(coregImage)
        if isfile(coregImage{i})
            ea_gencheckregpair(coregImage{i}, anchorImage, checkregFigure{i});
        end
    end
end

if ismember('brainshift', type)
    % Make dir
    ea_mkdir(fullfile(options.subj.subjDir, 'brainshift', 'checkreg'));

    % Get anchor Image path
    anchorImage = options.subj.brainshift.anat.anchor;

    % Generate checkreg figure for standard image
    if isfile(options.subj.brainshift.anat.moving)
        ea_gencheckregpair(options.subj.brainshift.anat.moving, anchorImage, options.subj.brainshift.checkreg.moving);
    end

    % Generate checkreg figure for brain shift corrected image
    if isfile(options.subj.brainshift.anat.scrf)
        ea_gencheckregpair(options.subj.brainshift.anat.scrf, anchorImage, options.subj.brainshift.checkreg.scrf);
    end
end

if ismember('norm', type)
    % Make dir
    ea_mkdir(fullfile(options.subj.subjDir, 'normalization', 'checkreg'));

    % Get template Image path
    templateImage = [ea_space, options.primarytemplate, '.nii'];

    % Remove CT from subj.norm.anat.postop struct
    if options.modality
        postop = rmfield(options.subj.norm.anat.postop, 'CT');
    end

    % Get paths of coregistered image
    normImage = [struct2cell(options.subj.norm.anat.preop); struct2cell(postop)];

    % Get paths of output figures
    checkregFigure = [struct2cell(options.subj.norm.checkreg.preop); struct2cell(options.subj.norm.checkreg.postop)];

    % Generate checkreg figures
    for i=1:length(normImage)
        if isfile(normImage{i})
            ea_gencheckregpair(normImage{i}, templateImage, checkregFigure{i});
        end
    end
end
