function options = ea_ensure_b0_coreg(options)
% EA_ENSURE_B0_COREG  Ensure that a coregistration transform between B0 and T1 exists.
% 
% This helper is intended to be called from the main pipeline (ea_autocoord)
% before any Lead-Connectome structural steps run. It is read-only with
% respect to options except for potential updates to options.coregmr.method
% and an in-memory LC-only correction of options.prefs.prenii_unnormalized
% via ea_lc_resolve_anat_anchor.
%
% Behaviour:
%   1. If a B0 image or preprocessed T1 cannot be found, the function
%      returns without doing anything (higher level code will error later).
%   2. It searches common locations for an existing transform B0 -> T1.
%   3. If none is found, it runs ea_coregimages once to create a
%      B0->T1 transform and stores it alongside the images.

directory = [options.root, options.patientname, filesep];

% Locate B0 image
b0Path = fullfile(directory, options.prefs.b0);
if ~isfile(b0Path)
    % Try BIDS-style location in preprocessing/dwi
    b0Dir = fullfile(directory, 'preprocessing', 'dwi');
    d = dir(fullfile(b0Dir, '*_b0.nii'));
    if isempty(d), d = dir(fullfile(b0Dir, '*b0*.nii')); end
    if isempty(d)
        fprintf('ea_ensure_b0_coreg: No B0 image found, skipping B0<->T1 coreg.\n');
        return;
    end
    b0Path = fullfile(d(1).folder, d(1).name);
end

% Locate anatomical anchor (LC-only resolver; ignores SPM/coreg intermediates)
[options, anatRel, anatName] = ea_lc_resolve_anat_anchor(options);
if isempty(anatRel)
    fprintf('ea_ensure_b0_coreg: No valid anatomical anchor found, skipping.\n');
    return;
end
anatPath = fullfile(directory, anatRel);
if ~isfile(anatPath)
    fprintf('ea_ensure_b0_coreg: Anatomical file missing (%s), skipping.\n', anatPath);
    return;
end

% Derive short names for pattern matching
[~, b0Name] = fileparts(b0Path);
b0Name  = regexprep(b0Name, '\.nii(\.gz)?$', '');
anatName = regexprep(anatName, '\.nii(\.gz)?$', '');

% Search for existing transform B0->T1 (include preprocessing/dwi where SPM often writes)
searchDirs = {
    fullfile(directory, 'coregistration', 'transformations')
    fullfile(directory, 'preprocessing', 'anat')
    fullfile(directory, 'preprocessing', 'dwi')
    directory
};
fwdPrefix = [b0Name, '2', anatName, '_'];
fwdAnts = [b0Name, '2', anatName, 'Composite'];
for iDir = 1:numel(searchDirs)
    cdir = searchDirs{iDir};
    if ~isfolder(cdir), continue; end
    mfiles = [dir(fullfile(cdir, '*.mat')); dir(fullfile(cdir, '*.h5')); dir(fullfile(cdir, '*Composite*.nii.gz'))];
    for k = 1:numel(mfiles)
        nm = mfiles(k).name;
        if startsWith(nm, fwdAnts)
            fprintf('ea_ensure_b0_coreg: Found existing B0->T1 transform: %s (in %s)\n', nm, cdir);
            return;
        end
        if startsWith(nm, fwdPrefix) && ~contains(nm, '_seg8') && ...
                (contains(nm, '_spm') || contains(lower(nm), 'ants') || contains(nm, 'Composite') || endsWith(nm, '.h5'))
            fprintf('ea_ensure_b0_coreg: Found existing B0->T1 transform: %s (in %s)\n', nm, cdir);
            return;
        end
    end
end

% No existing transform found -> run coregistration once
fprintf('ea_ensure_b0_coreg: No B0->T1 transform found. Running coregistration now...\n');

% Build output filename in preprocessing/anat
methodTok = regexp(options.coregmr.method, '^[^\s\(]+', 'match', 'once');
if isempty(methodTok), methodTok = 'spm'; end
outDir  = fullfile(directory, 'preprocessing', 'anat');
if ~isfolder(outDir), outDir = directory; end
ofile   = fullfile(outDir, [b0Name, '2', anatName, '.nii']);

try
    affinefile = ea_coregimages(options, b0Path, anatPath, ofile, {}, 1, [], 1);
    if ~isempty(affinefile)
        fprintf('ea_ensure_b0_coreg: Created B0->T1 transform: %s\n', affinefile{1});
    else
        fprintf('ea_ensure_b0_coreg: ea_coregimages did not return a transform file.\n');
    end
catch ME
    warning('ea_ensure_b0_coreg: Coregistration B0->T1 failed: %s', ME.message);
end
