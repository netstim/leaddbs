function ea_bids_tidy_subject_root(subjectDir)
% Move Lead-DBS/BIDS subject root files into the correct subfolders.
% Run once on a subject folder that had files written to root before
% the pipeline was updated to use BIDS folders.
%
% Example: ea_bids_tidy_subject_root('/path/to/derivatives/leaddbs/sub-XXX')

if nargin < 1 || isempty(subjectDir)
    subjectDir = uigetdir(pwd, 'Select subject directory (e.g. sub-XXX)');
    if isequal(subjectDir, 0)
        return;
    end
end
subjectDir = strtrim(subjectDir);
if subjectDir(end) == filesep
    subjectDir = subjectDir(1:end-1);
end
if ~isfolder(subjectDir)
    error('Not a directory: %s', subjectDir);
end

moved = 0;

% Normalization params -> normalization/transformations/
normDir = fullfile(subjectDir, 'normalization', 'transformations');
for f = {'y_ea_inv_normparams.nii', 'y_ea_normparams.nii'}
    src = fullfile(subjectDir, f{1});
    if isfile(src)
        if ~isfolder(normDir), mkdir(normDir); end
        movefile(src, fullfile(normDir, f{1}));
        fprintf('Moved %s -> normalization/transformations/\n', f{1});
        moved = moved + 1;
    end
end

% Tracking mask -> preprocessing/dwi/
dwiDir = fullfile(subjectDir, 'preprocessing', 'dwi');
src = fullfile(subjectDir, 'ttrackingmask.nii');
if isfile(src)
    if ~isfolder(dwiDir), mkdir(dwiDir); end
    movefile(src, fullfile(dwiDir, 'ttrackingmask.nii'));
    fprintf('Moved ttrackingmask.nii -> preprocessing/dwi/\n');
    moved = moved + 1;
end

% TR.mat -> prefs/
prefsDir = fullfile(subjectDir, 'prefs');
src = fullfile(subjectDir, 'TR.mat');
if isfile(src)
    if ~isfolder(prefsDir), mkdir(prefsDir); end
    movefile(src, fullfile(prefsDir, 'TR.mat'));
    fprintf('Moved TR.mat -> prefs/\n');
    moved = moved + 1;
end

% FA in anat: rename old _desc-FA_dwi.nii to _dwi_fa.nii
coregAnatDir = fullfile(subjectDir, 'coregistration', 'anat');
oldFa = dir(fullfile(coregAnatDir, '*_space-anchorNative_desc-FA_dwi.nii'));
for i = 1:numel(oldFa)
    src = fullfile(oldFa(i).folder, oldFa(i).name);
    newName = strrep(oldFa(i).name, '_desc-FA_dwi.nii', '_dwi_fa.nii');
    if ~strcmp(oldFa(i).name, newName)
        movefile(src, fullfile(oldFa(i).folder, newName));
        fprintf('Renamed %s -> %s\n', oldFa(i).name, newName);
        moved = moved + 1;
    end
end

% Coregistration method log -> coregistration/log/
coregLogDir = fullfile(subjectDir, 'coregistration', 'log');
src = fullfile(subjectDir, 'ea_coregmrmethod_applied.mat');
if isfile(src)
    if ~isfolder(coregLogDir), mkdir(coregLogDir); end
    movefile(src, fullfile(coregLogDir, 'ea_coregmrmethod_applied.mat'));
    fprintf('Moved ea_coregmrmethod_applied.mat -> coregistration/log/\n');
    moved = moved + 1;
end

% Resliced/SPM outputs (r*.*.nii) -> preprocessing/anat/ (legacy root cleanup)
anatDir = fullfile(subjectDir, 'preprocessing', 'anat');
d = dir(fullfile(subjectDir, 'r*.nii'));
if isempty(d), d = struct('folder', {}, 'name', {}); end
for i = 1:numel(d)
    src = fullfile(d(i).folder, d(i).name);
    if isfile(src)
        if ~isfolder(anatDir), mkdir(anatDir); end
        movefile(src, fullfile(anatDir, d(i).name));
        fprintf('Moved %s -> preprocessing/anat/\n', d(i).name);
        moved = moved + 1;
    end
end

if moved == 0
    fprintf('No root-level files to move in %s\n', subjectDir);
else
    fprintf('Done. Moved %d item(s).\n', moved);
end
