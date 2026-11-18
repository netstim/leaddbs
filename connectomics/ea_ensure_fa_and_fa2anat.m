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

% If there is already a FA coregistration to anat, reuse it unless it is
% older than the coregistered B0 (i.e. the B0 transform changed since).
if isfile(fa2anatPath)
    if isBIDS
        options.prefs.fa2anat = fa2anatRel;
    end
    % Always add the existing FA to the options, so ea_gencheckregfigs and
    % ea_checkreg (called next in ea_autocoord) include it. This matters on
    % a first run: ea_getptopts ran before the coregistration, when the
    % FA in anchor space (created together with the B0) did not exist yet.
    options = add_fa_to_options(options, fa2anatPath, faPath);
    if ~fa_is_stale(fa2anatPath, options)
        return;
    end
    fprintf('ea_ensure_fa_and_fa2anat: FA in anchor space is older than the coregistered B0. Rebuilding it.\n');
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
[transform, b0Method] = find_b0_t1_forward_transform(directory, options);
if isempty(transform)
    warning(['ea_ensure_fa_and_fa2anat: No saved B0->anchor transform for the method used for the B0 ("%s"). ', ...
             'The FA cannot be rebuilt on its own. Please rerun the B0 coregistration ', ...
             '(the FA is carried along with it).'], b0Method);
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
    % ea_checkreg (called next in ea_autocoord) show FA on this same run.
    options = add_fa_to_options(options, fa2anatPath, faPath);
catch ME
    warning('ea_ensure_fa_and_fa2anat: Failed to apply transform to FA: %s', ME.message);
end


function options = add_fa_to_options(options, fa2anatPath, faPath)
% Add the FA (anchor space) to the coreg/preproc/checkreg structs, the same
% way ea_getptopts does, so it gets a checkreg figure and is reviewed.
if ~isfile(fa2anatPath) || ...
        ~isfield(options, 'subj') || ...
        ~isfield(options.subj, 'coreg') || ...
        ~isfield(options.subj.coreg, 'anat') || ...
        ~isfield(options.subj.coreg.anat, 'preop') || ...
        ~isfield(options.subj, 'preproc') || ...
        ~isfield(options.subj.preproc, 'anat') || ...
        ~isfield(options.subj.preproc.anat, 'preop')
    return;
end

% Already the DWI-derived FA: nothing to do. Any other 'fa' entry (e.g. an FA
% imported as anatomical image, whose anchor-space file never exists) is
% replaced by the DWI-derived FA, which carries the B0 transform.
if isfield(options.subj.coreg.anat.preop, 'fa') && ...
        strcmp(options.subj.coreg.anat.preop.fa, fa2anatPath)
    return;
end

options.subj.coreg.anat.preop.fa   = fa2anatPath;
options.subj.preproc.anat.preop.fa = faPath;

if isfield(options.subj, 'coregDir') && ...
        isfield(options.subj.coreg, 'checkreg') && ...
        isfield(options.subj.coreg.checkreg, 'preop')
    [~, faCoregName] = fileparts(fa2anatPath);
    options.subj.coreg.checkreg.preop.fa = fullfile( ...
        options.subj.coregDir, 'checkreg', [faCoregName, '.png']);
end

function [hit, method] = find_b0_t1_forward_transform(directory, options) %#ok<INUSD>
% Return the B0->anchor forward transform written by the B0 coregistration,
% or '' if the method used for the B0 does not save a reusable transform.
%
% The method is taken from the coregistration log (method.b0), i.e. the
% method that actually produced the current B0 — not the method currently
% selected in the GUI. Transform files left over from other methods are
% never used.

hit = '';

% 1. Method actually used for the B0
method = '';
if isfield(options, 'subj') && isfield(options.subj, 'coreg') && ...
        isfield(options.subj.coreg, 'log') && isfile(options.subj.coreg.log.method)
    try
        json = loadjson(options.subj.coreg.log.method);
        if isfield(json, 'method') && isfield(json.method, 'b0')
            method = json.method.b0;
        end
    catch
    end
end
if isempty(method) && isfield(options, 'coregmr') && isfield(options.coregmr, 'method')
    method = options.coregmr.method;
end

% 2. Suffix written by ea_coregpreopmr for this method. Methods without a
% saved transform (FLIRT BBR, hybrids, ANTs nonlinear) return ''.
switch lower(method)
    case lower({'ANTs (Avants 2008)', 'ANTs'})
        suffix = 'ants.mat';
    case lower({'FLIRT (Jenkinson 2001 & 2002)', 'FLIRT'})
        suffix = 'flirt.mat';
    case lower({'SPM (Friston 2007)', 'SPM'})
        suffix = 'spm.mat';
    case lower({'BRAINSFit (Johnson 2007)', 'BRAINSFit'})
        suffix = 'brainsfit.mat';
    otherwise
        return;
end

if ~isfield(options.subj.coreg, 'transform') || ~isfield(options.subj.coreg.transform, 'b0')
    return;
end
candidate = [options.subj.coreg.transform.b0.forwardBaseName, suffix];
if isfile(candidate)
    hit = candidate;
end


function stale = fa_is_stale(fa2anatPath, options)
% FA in anchor space is stale when it is older than the coregistered B0.
stale = false;
if ~isfield(options, 'subj') || ~isfield(options.subj, 'coreg') || ...
        ~isfield(options.subj.coreg.anat, 'preop') || ...
        ~isfield(options.subj.coreg.anat.preop, 'b0') || ...
        ~isfile(options.subj.coreg.anat.preop.b0)
    return;
end
faInfo = dir(fa2anatPath);
b0Info = dir(options.subj.coreg.anat.preop.b0);
% 2 s tolerance for file system time stamp resolution
stale = faInfo.datenum < b0Info.datenum - 2/86400;
