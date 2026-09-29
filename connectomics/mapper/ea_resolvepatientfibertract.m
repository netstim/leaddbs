function [fiberFile, candidates] = ea_resolvepatientfibertract(patientDir, useNativeSeed, prefs)
% Resolve a patient-specific fiber tract file across BIDS and legacy layouts.

if nargin < 2 || isempty(useNativeSeed)
    useNativeSeed = false;
end
if nargin < 3 || isempty(prefs)
    prefs = struct;
end
patientDir = char(patientDir);

if useNativeSeed
    bidsFile = 'FTR_anat.mat';
    if isfield(prefs, 'FTR_unnormalized')
        prefFile = prefs.FTR_unnormalized;
    else
        prefFile = 'FTR.mat';
    end
    [rawPrefDir, rawPrefName, rawPrefExt] = fileparts(prefFile);
    isDefaultPreference = strcmp([rawPrefName, rawPrefExt], 'FTR.mat') && ...
        (isempty(rawPrefDir) || strcmp(rawPrefDir, fullfile('connectomics', 'dMRI')));
    [prefDir, prefName, prefExt] = fileparts(prefFile);
    if ~endsWith(prefName, '_anat')
        prefName = [prefName, '_anat'];
    end
    prefFile = fullfile(prefDir, [prefName, prefExt]);
else
    bidsFile = 'FTR_normalized.mat';
    if isfield(prefs, 'FTR_normalized')
        prefFile = prefs.FTR_normalized;
    else
        prefFile = 'wFTR.mat';
    end
    isDefaultPreference = any(strcmp(char(prefFile), ...
        {'wFTR.mat', fullfile('connectomics', 'dMRI', 'FTR_normalized.mat')}));
end

prefFile = char(prefFile);
[~, prefName, prefExt] = fileparts(prefFile);
prefBaseName = [prefName, prefExt];

if startsWith(prefFile, filesep) || ...
        ~isempty(regexp(prefFile, '^[A-Za-z]:[\\/]', 'once')) || ...
        startsWith(prefFile, '\\')
    preferredFile = prefFile;
else
    preferredFile = fullfile(patientDir, prefFile);
end

bidsFile = fullfile(patientDir, 'connectomics', 'dMRI', bidsFile);
fallbacks = {
    preferredFile
    fullfile(patientDir, 'connectomics', 'dMRI', prefBaseName)
    fullfile(patientDir, 'connectomes', 'dMRI', prefBaseName)
    fullfile(patientDir, prefBaseName)
    };

% Conventional preferences still point to legacy names, so prefer the BIDS
% path for those. An explicit custom preference should retain precedence.
if isDefaultPreference
    candidates = [{bidsFile}; fallbacks];
else
    candidates = [fallbacks; {bidsFile}];
end
candidates = unique(candidates, 'stable');

fiberFile = '';
existingFile = find(isfile(candidates), 1);
if ~isempty(existingFile)
    fiberFile = candidates{existingFile};
end
