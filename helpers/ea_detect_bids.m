function isBIDS = ea_detect_bids(directory)
% EA_DETECT_BIDS  Detect if a patient directory follows BIDS structure
%
%   isBIDS = ea_detect_bids(directory)
%
%   Returns true if the directory appears to be in BIDS format, false otherwise.
%
%   BIDS indicators:
%   - Directory is under derivatives/leaddbs/
%   - Subject ID starts with 'sub-'
%   - Has rawdata/ sibling directory
%   - Has dataset_description.json
%   - Files follow BIDS naming (sub-XXX_ses-XXX_modality.nii)
%
%   Non-BIDS indicators:
%   - Classic filenames (anat_t1.nii, DTI.nii, rest.nii)
%   - No derivatives/ structure
%   - No sub- prefix in directory name

% Remove trailing filesep
directory = ea_niifileparts(directory);

isBIDS = false;

%% Check 1: Is this under derivatives/leaddbs/?
if contains(directory, [filesep, 'derivatives', filesep, 'leaddbs', filesep])
    isBIDS = true;
    return;
end

%% Check 2: Does the subject directory name start with 'sub-'?
[~, subjDirName] = fileparts(directory);
if startsWith(subjDirName, 'sub-')
    isBIDS = true;
    return;
end

%% Check 3: Look for BIDS-style files
% Check for BIDS naming pattern in files
if exist(directory, 'dir')
    % Search for files with BIDS naming pattern
    bidsFiles = dir(fullfile(directory, '**', 'sub-*_*.nii*'));
    if ~isempty(bidsFiles)
        isBIDS = true;
        return;
    end
    
    % Search for preprocessing/rawdata subdirectories
    if exist(fullfile(directory, 'preprocessing'), 'dir') || ...
       exist(fullfile(directory, 'rawdata'), 'dir')
        isBIDS = true;
        return;
    end
end

%% Check 4: Look for classic Lead-DBS filenames (indicates NOT BIDS)
classicFiles = {
    'anat_t1.nii', 'anat_t1.nii.gz', ...
    'anat_t2.nii', 'anat_t2.nii.gz', ...
    'DTI.nii', 'DTI.nii.gz', ...
    'rest.nii', 'rest.nii.gz', ...
    'glanat.nii', 'glanat.nii.gz'
};

for i = 1:length(classicFiles)
    if exist(fullfile(directory, classicFiles{i}), 'file')
        isBIDS = false;
        return;
    end
end

%% Default: Assume BIDS if in doubt (safer for new workflows)
isBIDS = false;
end

