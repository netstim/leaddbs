function [tra, presentfiles] = ea_assignpretra(options)
% Stub function for BIDS compatibility
% Returns anatomical files present in the patient directory

directory = [options.root, options.patientname, filesep];

% Check for BIDS preprocessing/anat files
preprocAnatDir = fullfile(directory, 'preprocessing', 'anat');
if exist(preprocAnatDir, 'dir')
    % BIDS structure: look for T1w, T2w in preprocessing/anat
    anatFiles = dir(fullfile(preprocAnatDir, '*_T1w.nii'));
    if isempty(anatFiles)
        anatFiles = dir(fullfile(preprocAnatDir, '*_T2w.nii'));
    end
    if ~isempty(anatFiles)
        presentfiles = {fullfile('preprocessing', 'anat', anatFiles(1).name)};
        tra = presentfiles;
        return;
    end
end

% Fallback: classic file structure
tra = {};
presentfiles = {};

% Check for classic anatomical files
if exist(fullfile(directory, 'anat_t1.nii'), 'file')
    tra{end+1} = 'anat_t1.nii';
    presentfiles{end+1} = 'anat_t1.nii';
elseif exist(fullfile(directory, 'anat_t2.nii'), 'file')
    tra{end+1} = 'anat_t2.nii';
    presentfiles{end+1} = 'anat_t2.nii';
end

% If still empty, use prenii_unnormalized from options
if isempty(presentfiles) && isfield(options, 'prefs') && isfield(options.prefs, 'prenii_unnormalized')
    presentfiles = {options.prefs.prenii_unnormalized};
    tra = presentfiles;
end

