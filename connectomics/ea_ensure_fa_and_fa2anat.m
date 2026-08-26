function options = ea_ensure_fa_and_fa2anat(options)
% Ensure FA map exists and FA coregistered to anatomical is in coregistration/anat.
%
% When running Lead Connectome (structural), this helper:
%  1. Creates FA from DWI if not present (preprocessing/dwi/*_fa.nii).
%  2. Coregisters FA to the anatomical (T1) and writes the result to
%     coregistration/anat (BIDS) or subject root fa2anat.nii (legacy).
%
% Called from ea_autocoord when any structural connectome option is enabled.

directory = [options.root, options.patientname, filesep];

% Need DWI data (options.prefs.dti set by ea_prepare_dti_bids or legacy)
if ~isfield(options.prefs, 'dti') || isempty(options.prefs.dti)
    return;
end
dtiPath = fullfile(directory, options.prefs.dti);
if ~isfile(dtiPath)
    return;
end

% 1) Create FA if missing
faPath = fullfile(directory, options.prefs.fa);
if ~isfile(faPath)
    fprintf('\nCreating FA map from DWI...\n');
    try
        ea_isolate_fa(options);
        fprintf('FA saved: %s\n', options.prefs.fa);
    catch ME
        warning('ea_ensure_fa_and_fa2anat: Could not create FA: %s', ME.message);
        return;
    end
end

% 2) FA-in-anat: output path
isBIDS = contains(directory, 'derivatives') || contains(directory, 'leaddbs');
if isBIDS
    coregAnatDir = fullfile(directory, 'coregistration', 'anat');
    if ~isfolder(coregAnatDir)
        mkdir(coregAnatDir);
    end
    % BIDS-style name: sub-XXX_space-anchorNative_dwi_fa.nii
    fa2anatName = [options.patientname, '_space-anchorNative_dwi_fa.nii'];
    fa2anatPath = fullfile(coregAnatDir, fa2anatName);
    fa2anatRel  = fullfile('coregistration', 'anat', fa2anatName);
else
    fa2anatPath = fullfile(directory, options.prefs.fa2anat);
    fa2anatRel  = options.prefs.fa2anat;
end

if isfile(fa2anatPath)
    if isBIDS
        options.prefs.fa2anat = fa2anatRel;
    end
    return;
end

% Anatomical reference (LC-only resolver; ignores SPM/coreg intermediates)
[options, anatRel] = ea_lc_resolve_anat_anchor(options);
if isempty(anatRel)
    warning('ea_ensure_fa_and_fa2anat: Anatomical reference not found. Skipping FA->anat coregistration.');
    return;
end
anatPath = fullfile(directory, anatRel);
if ~isfile(anatPath)
    warning('ea_ensure_fa_and_fa2anat: Anatomical reference not found (%s). Skipping FA->anat coregistration.', anatPath);
    return;
end

% Run FA -> anat coregistration; write to fa2anat path
fprintf('\nCoregistering FA to anatomical and saving to coregistration/anat...\n');
try
    ea_coregimages(options, faPath, anatPath, fa2anatPath, {}, 1, [], 1);
    if isBIDS
        options.prefs.fa2anat = fa2anatRel;
    end
    fprintf('FA-in-anat saved: %s\n', fa2anatRel);
catch ME
    warning('ea_ensure_fa_and_fa2anat: FA->anat coregistration failed: %s', ME.message);
end
