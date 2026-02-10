function options = ea_ensure_b0_coreg(options)
% EA_ENSURE_B0_COREG  Ensure that a coregistration transform between B0 and T1 exists.
% 
% This helper is intended to be called from the main pipeline (ea_autocoord)
% before any Lead-Connectome structural steps run. It is read-only with
% respect to options except for potential updates to options.coregmr.method.
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

% Locate preprocessed T1
anatDir = fullfile(directory, 'preprocessing', 'anat');
if ~isfolder(anatDir)
    fprintf('ea_ensure_b0_coreg: preprocessing/anat directory not found, skipping.\n');
    return;
end
cand = dir(fullfile(anatDir, '*desc-preproc*_T1w.nii'));
if isempty(cand)
    cand = dir(fullfile(anatDir, '*acq-iso_T1w.nii'));
end
if isempty(cand)
    cand = dir(fullfile(anatDir, '*T1w.nii'));
end
if isempty(cand)
    fprintf('ea_ensure_b0_coreg: No preprocessed T1 found in %s, skipping.\n', anatDir);
    return;
end
anatPath = fullfile(cand(1).folder, cand(1).name);

% Derive short names for pattern matching
[~, b0Name] = fileparts(b0Path);
[~, anatName] = fileparts(anatPath);
b0Name  = regexprep(b0Name, '\.nii(\.gz)?$', '');
anatName = regexprep(anatName, '\.nii(\.gz)?$', '');

% Search for existing transform B0->T1
searchDirs = {
    fullfile(directory, 'coregistration', 'transformations')
    fullfile(directory, 'preprocessing', 'anat')
    directory
};
pattern = [b0Name, '2', anatName]; % B0 -> T1
for iDir = 1:numel(searchDirs)
    cdir = searchDirs{iDir};
    if ~isfolder(cdir), continue; end
    mfiles = dir(fullfile(cdir, '*.mat'));
    for k = 1:numel(mfiles)
        if contains(mfiles(k).name, pattern)
            fprintf('ea_ensure_b0_coreg: Found existing B0->T1 transform: %s (in %s)\n', ...
                mfiles(k).name, cdir);
            return;
        end
    end
end

% No existing transform found -> run coregistration once
fprintf('ea_ensure_b0_coreg: No B0->T1 transform found. Running coregistration now...\n');

% Build output filename in preprocessing/anat
outName = sprintf('%s2%s_%s.mat', b0Name, anatName, lower(regexp(options.coregmr.method, '^[^\s\(]+', 'match', 'once')));
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

