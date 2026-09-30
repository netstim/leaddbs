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
    ea_mkdir(coregAnatDir);
    % BIDS-style name: sub-XXX_ses-preop_space-anchorNative_fa.nii
    fa2anatName = ['sub-', options.subj.subjId, '_ses-preop_space-anchorNative_fa.nii'];
    fa2anatPath = fullfile(coregAnatDir, fa2anatName);
    fa2anatRel  = fullfile('coregistration', 'anat', fa2anatName);
else
    fa2anatPath = fullfile(directory, options.prefs.fa2anat);
    fa2anatRel  = options.prefs.fa2anat;
end

% If there is already a FA coregistration to anat
if isfile(fa2anatPath)
    if isBIDS
        options.prefs.fa2anat = fa2anatRel;
    end
    return;
end

% Use the same anchor normalization will use; prenii_unnormalized can
% point to a different grid and cause a mismatch later.
anatPath = '';
if isfield(options.subj.coreg.anat.preop, options.subj.AnchorModality)
    candidate = options.subj.coreg.anat.preop.(options.subj.AnchorModality);
    if isfile(candidate)
        anatPath = candidate;
    end
end
if isempty(anatPath)
    % Legacy fallback
    candidate = fullfile(directory, options.prefs.prenii_unnormalized);
    if isfile(candidate)
        anatPath = candidate;
    else
        for subdir = {'preprocessing/anat', 'coregistration/anat'}
            d = dir(fullfile(directory, subdir{1}, '*T1w.nii'));
            if isempty(d), d = dir(fullfile(directory, subdir{1}, '*T2w.nii')); end
            if ~isempty(d)
                anatPath = fullfile(d(1).folder, d(1).name);
                break;
            end
        end
    end
end
if isempty(anatPath) || ~isfile(anatPath)
    warning('ea_ensure_fa_and_fa2anat: Anatomical reference not found. Skipping FA->anat coregistration.');
    return;
end
fprintf('ea_ensure_fa_and_fa2anat: Using anchor reference: %s\n', anatPath);

% Find the B0->T1 forward transform 
transform = find_b0_t1_forward_transform(directory, options);
if isempty(transform)
    warning(['ea_ensure_fa_and_fa2anat: B0->T1 forward transform not found. ', ...
             'Ensure B0 coregistration ran successfully.']);
    return;
end
fprintf('ea_ensure_fa_and_fa2anat: Using B0->T1 transform: %s\n', transform);

% Apply transform to FA
fprintf('ea_ensure_fa_and_fa2anat: Applying B0->T1 transform to FA...\n');
try
    ea_apply_coregistration(anatPath, faPath, fa2anatPath, transform);

    if isBIDS
        options.prefs.fa2anat = fa2anatRel;
    end

    % Inject FA into the options structs so that ea_gencheckregfigs and
    % ea_checkreg (called next in ea_autocoord) show FA on this same run,
    % without waiting for a second call to ea_getptopts.
    if isfield(options, 'subj') && ...
            isfield(options.subj, 'coreg') && ...
            isfield(options.subj.coreg, 'anat') && ...
            isfield(options.subj.coreg.anat, 'preop') && ...
            ~isfield(options.subj.coreg.anat.preop, 'fa') && ...
            isfield(options.subj, 'preproc') && ...
            isfield(options.subj.preproc, 'anat') && ...
            isfield(options.subj.preproc.anat, 'preop')

        options.subj.coreg.anat.preop.fa   = fa2anatPath;
        options.subj.preproc.anat.preop.fa = faPath;

        if isfield(options.subj, 'coregDir') && ...
                isfield(options.subj.coreg, 'checkreg') && ...
                isfield(options.subj.coreg.checkreg, 'preop')
            [~, faCoregName] = fileparts(fa2anatPath);
            options.subj.coreg.checkreg.preop.fa = fullfile( ...
                options.subj.coregDir, 'checkreg', [faCoregName, '.png']);
        end
    end
catch ME
    warning('ea_ensure_fa_and_fa2anat: Failed to apply transform to FA: %s', ME.message);
end

function hit = find_b0_t1_forward_transform(directory, options)
% Return full path to a B0->T1 forward transform file, or '' if none found.
%
% Strategy:
%  1. Use options.subj.coreg.transform.b0.forwardBaseName directly
%     (most reliable — already populated by ea_coregpreopmr).
%  2. Fall back to a scored file-system search in coregistration/transformations/.

hit = '';

% Extract method string
if isstruct(options) && isfield(options, 'coregmr') && isfield(options.coregmr, 'method')
    method = options.coregmr.method;
elseif ischar(options)
    method = options;
else
    method = '';
end
methodHint = lower(regexp(method, '^[^\s\(]+', 'match', 'once'));

% ── 1. Struct-based lookup (preferred) ───────────────────────────────────
if isstruct(options) && isfield(options, 'subj') && ...
        isfield(options.subj, 'coreg') && ...
        isfield(options.subj.coreg, 'transform') && ...
        isfield(options.subj.coreg.transform, 'b0')

    base = options.subj.coreg.transform.b0.forwardBaseName;

    % Map method string to the transform file suffix saved by ea_coregpreopmr
    switch methodHint
        case 'spm'
            suffixes = {'spm.mat'};
        case 'ants'
            suffixes = {'ants.mat'};   % ITK affine – correct input for antsApplyTransforms
        case {'flirt', 'flirtbbr', 'bbr', 'fsl'}
            suffixes = {'flirt.mat'};
        case 'brainsfit'
            suffixes = {'brainsfit.mat'};
        otherwise
            suffixes = {'spm.mat', 'ants.mat', 'flirt.mat', 'brainsfit.mat'};
    end

    for k = 1:numel(suffixes)
        candidate = [base, suffixes{k}];
        if isfile(candidate)
            hit = candidate;
            return;
        end
    end
end

% ── 2. File-system search fallback ───────────────────────────────────────
searchDir = fullfile(directory, 'coregistration', 'transformations');
if ~isfolder(searchDir), return; end

exts  = {'*.mat', '*.h5', '*.txt'};
cands = {};
for e = 1:numel(exts)
    d = dir(fullfile(searchDir, '**', exts{e}));
    for k = 1:numel(d)
        cands{end+1} = fullfile(d(k).folder, d(k).name); %#ok<AGROW>
    end
end
if isempty(cands), return; end

bestScore = -Inf;
for i = 1:numel(cands)
    [~, name, ext] = fileparts(cands{i});
    fname = lower([name, ext]);

    % Must reference both B0/DWI and T1/anat side
    hasB0   = contains(fname, 'b0') || contains(fname, 'dwi');
    hasAnat = contains(fname, 't1') || contains(fname, 'anat') || ...
              contains(fname, 'anchor') || contains(fname, 'native');
    if ~hasB0 || ~hasAnat, continue; end

    % Exclude inverse transforms
    isInverse = startsWith(fname, 'anat') || startsWith(fname, 't1') || ...
                contains(fname, 'from-anchor') || contains(fname, 'from-t1');
    if isInverse, continue; end

    % Exclude the ants44 / spm44 / flirt44 convenience copies — those are
    % MATLAB-format 4x4 matrices, not suitable inputs for the apply functions.
    if regexp(fname, '\d+\.mat$'), continue; end

    score = 4; % baseline

    if ~isempty(methodHint) && contains(fname, methodHint)
        score = score + 1;
    end

    % For MATLAB-native .mat files (SPM), verify they contain a 4x4 matrix.
    % ANTs ITK .mat files are binary and cannot be loaded by MATLAB — don't
    % penalise them; their name already identifies them.
    if strcmp(ext, '.mat') && ~contains(fname, 'ants')
        try
            S = load(cands{i});
            has4x4 = any(structfun(@(v) isnumeric(v) && isequal(size(v), [4 4]), S));
            if has4x4
                score = score + 1;
            else
                score = score - 2;
            end
        catch
            score = score - 3;
        end
    end

    if score > bestScore
        bestScore = score;
        hit = cands{i};
    end
end

if bestScore < 4
    hit = '';
end
