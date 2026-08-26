function [options, anatRel, anatName] = ea_lc_resolve_anat_anchor(options)
% EA_LC_RESOLVE_ANAT_ANCHOR  Resolve anatomical anchor for Lead-Connectome only.
%
% Finds the real T1/T2 anatomical image for DWI<->anat coregistration and
% fiber normalization, ignoring SPM/coreg intermediate products that often
% pollute preprocessing/anat (c1*/w*/b02*/dwi* nested names).
%
% Preference order:
%   1. coregistration/anat space-anchorNative preproc T1/T2
%   2. preprocessing/anat desc-preproc T1/T2 (filtered)
%   3. existing prefs.prenii_unnormalized if it passes the filter
% Among valid candidates, prefer one that already has a B0->anat transform
% on disk (avoids unnecessary re-coreg when preprocessing anat was used).
%
% Does NOT modify ea_getptopts or any global prefs files. Only updates the
% in-memory options.prefs.prenii_unnormalized for the current LC call.
%
% Outputs:
%   options  - options with prenii_unnormalized corrected when a better
%              anchor was found
%   anatRel  - path relative to subject directory (or '' if unresolved)
%   anatName - basename without .nii/.nii.gz (or '' if unresolved)

directory = [options.root, options.patientname, filesep];
anatRel = '';
anatName = '';

candidates = {}; % full paths, preference order

% --- Priority 1: coregistration/anat anchorNative preproc ---
coregAnatDir = fullfile(directory, 'coregistration', 'anat');
candidates = [candidates, collect_valid(coregAnatDir, {
    '*space-anchorNative*desc-preproc*T1w.nii'
    '*space-anchorNative*T1w.nii'
    '*space-anchorNative*desc-preproc*T2w.nii'
    '*space-anchorNative*T2w.nii'
})];

% --- Priority 2: preprocessing/anat desc-preproc T1/T2 ---
preprocAnatDir = fullfile(directory, 'preprocessing', 'anat');
candidates = [candidates, collect_valid(preprocAnatDir, {
    '*desc-preproc*_T1w.nii'
    '*desc-preproc*_T2w.nii'
    '*_T1w.nii'
    '*_T2w.nii'
})];

% --- Priority 3: existing prefs.prenii_unnormalized if it passes filter ---
if isfield(options, 'prefs') && isfield(options.prefs, 'prenii_unnormalized') ...
        && ~isempty(options.prefs.prenii_unnormalized)
    prefPath = fullfile(directory, options.prefs.prenii_unnormalized);
    if isfile(prefPath)
        [~, prefBase, prefExt] = fileparts(options.prefs.prenii_unnormalized);
        if strcmpi(prefExt, '.gz')
            [~, prefBase] = fileparts(prefBase);
        end
        if is_valid_anat_basename(prefBase)
            candidates{end+1} = prefPath; %#ok<AGROW>
        end
    end
end

% Deduplicate while preserving order
candidates = unique_preserve(candidates);

if isempty(candidates)
    fprintf('ea_lc_resolve_anat_anchor: No valid anatomical anchor found.\n');
    return;
end

% Prefer candidate that already has a B0->anat transform (if b0 known)
cand = candidates{1};
b0Name = '';
if isfield(options, 'prefs') && isfield(options.prefs, 'b0') && ~isempty(options.prefs.b0)
    b0Path = fullfile(directory, options.prefs.b0);
    if ~isfile(b0Path)
        b0Dir = fullfile(directory, 'preprocessing', 'dwi');
        d = dir(fullfile(b0Dir, '*_b0.nii'));
        if isempty(d), d = dir(fullfile(b0Dir, '*b0*.nii')); end
        if ~isempty(d)
            b0Path = fullfile(d(1).folder, d(1).name);
        end
    end
    if isfile(b0Path)
        [~, b0Name] = fileparts(b0Path);
        b0Name = regexprep(b0Name, '\.nii(\.gz)?$', '');
    end
end

if ~isempty(b0Name)
    for i = 1:numel(candidates)
        [~, thisAnat] = fileparts(candidates{i});
        thisAnat = regexprep(thisAnat, '\.nii(\.gz)?$', '');
        if has_b0_to_anat_transform(directory, b0Name, thisAnat)
            cand = candidates{i};
            break;
        end
    end
end

% Relative path from subject root
if startsWith(cand, directory)
    anatRel = cand(length(directory)+1:end);
else
    anatRel = cand;
end
if filesep ~= '/'
    anatRel = strrep(strrep(anatRel, '\', '/'), '/', filesep);
else
    anatRel = strrep(anatRel, '\', '/');
end

[~, anatName] = fileparts(anatRel);
anatName = regexprep(anatName, '\.nii(\.gz)?$', '');

% In-memory correction for this LC run only
prev = '';
if isfield(options, 'prefs') && isfield(options.prefs, 'prenii_unnormalized')
    prev = options.prefs.prenii_unnormalized;
end
options.prefs.prenii_unnormalized = anatRel;

if ~strcmp(prev, anatRel)
    fprintf('ea_lc_resolve_anat_anchor: Using anatomical anchor: %s\n', anatRel);
    if ~isempty(prev)
        fprintf('  (prefs.prenii_unnormalized was: %s)\n', prev);
    end
end


function paths = collect_valid(searchDir, patterns)
paths = {};
if ~isfolder(searchDir)
    return;
end
for p = 1:numel(patterns)
    hits = dir(fullfile(searchDir, patterns{p}));
    for k = 1:numel(hits)
        if hits(k).isdir
            continue;
        end
        [~, base, ext] = fileparts(hits(k).name);
        if strcmpi(ext, '.gz')
            [~, base] = fileparts(base);
        end
        if is_valid_anat_basename(base)
            paths{end+1} = fullfile(hits(k).folder, hits(k).name); %#ok<AGROW>
        end
    end
end


function out = unique_preserve(in)
out = {};
for i = 1:numel(in)
    if ~any(strcmp(out, in{i}))
        out{end+1} = in{i}; %#ok<AGROW>
    end
end


function tf = has_b0_to_anat_transform(directory, b0Name, anatName)
tf = false;
% Forward coreg mats start with b0Name2anatName_ and carry a method token
% (e.g. _spm). Exclude SPM segmentation sidecars (*_seg8.mat).
fwdPrefix = [b0Name, '2', anatName, '_'];
fwdAnts = [b0Name, '2', anatName, 'Composite'];
searchDirs = {
    fullfile(directory, 'preprocessing', 'anat')
    fullfile(directory, 'preprocessing', 'dwi')
    fullfile(directory, 'coregistration', 'transformations')
    directory
};
for iDir = 1:numel(searchDirs)
    cdir = searchDirs{iDir};
    if ~isfolder(cdir), continue; end
    mfiles = [dir(fullfile(cdir, '*.mat')); dir(fullfile(cdir, '*.h5')); dir(fullfile(cdir, '*Composite*.nii.gz'))];
    for k = 1:numel(mfiles)
        nm = mfiles(k).name;
        if startsWith(nm, fwdAnts)
            tf = true;
            return;
        end
        if startsWith(nm, fwdPrefix) && ~contains(nm, '_seg8') && ...
                (contains(nm, '_spm') || contains(lower(nm), 'ants') || contains(nm, 'Composite') || endsWith(nm, '.h5'))
            tf = true;
            return;
        end
    end
end


function tf = is_valid_anat_basename(base)
% Reject SPM tissue classes, warped outputs, and moving2fixed intermediates.
tf = false;
if isempty(base)
    return;
end

% Prefix filters (c1/c2/c3, w, r, sr, mean)
if ~isempty(regexp(base, '^(c[1-3]|w|r|sr|mean)', 'once'))
    return;
end

% Coreg / DWI intermediate products (b02..., ...2c1..., ...2sub-...)
if contains(base, 'dwi', 'IgnoreCase', true)
    return;
end
if contains(base, 'b02')
    return;
end
if contains(base, '2c1') || contains(base, '2c2') || contains(base, '2c3')
    return;
end
if contains(base, '2sub-')
    return;
end

tf = true;
