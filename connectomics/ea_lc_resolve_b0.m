function [options, b0Rel, b0Path] = ea_lc_resolve_b0(options)
% Resolve the native B0 image used by tractography.
%
% BIDS subjects can also contain a resliced B0-to-anatomy image whose name
% includes "b02<anat>".  That image must never be used as the source grid
% for native FTR coordinates.

directory = fullfile(options.root, options.patientname);
b0Path = '';

% Prefer the canonical BIDS preprocessing product over a stale prefs value.
b0Files = [dir(fullfile(directory, 'preprocessing', 'dwi', '*_b0.nii')); ...
           dir(fullfile(directory, 'preprocessing', 'dwi', '*_b0.nii.gz'))];
if ~isempty(b0Files)
    names = {b0Files.name};
    keep = ~contains(names, 'b02') & ~startsWith(names, {'r', 'w'});
    b0Files = b0Files(keep);
end
if ~isempty(b0Files)
    [~, newest] = max([b0Files.datenum]);
    b0Path = fullfile(b0Files(newest).folder, b0Files(newest).name);
end

if isempty(b0Path) && isfield(options.prefs, 'b0') && ~isempty(options.prefs.b0)
    candidate = options.prefs.b0;
    if ~isfile(candidate)
        candidate = fullfile(directory, candidate);
    end
    [~, candidateName] = ea_niifileparts(candidate);
    if isfile(candidate) && ~contains(candidateName, 'b02')
        b0Path = candidate;
    end
end

if isempty(b0Path)
    b0Rel = '';
    return;
end

b0Rel = erase(b0Path, [directory, filesep]);
options.prefs.b0 = b0Rel;

