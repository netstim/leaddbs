function ea_dataset_import(source_dir, dest_dir, method, dicomimport)

% This function has two functionalities:
%
% (1) import an already BIDS compliant dataset (with rawdata available):
% screens the rawdata folder for appropriate pre-
% and postop files. If more than one file is found, a small window is
% triggerd and the user has to select the file they want lead-dbs to use
% A typical call would be ea_dataset_import(path_to_root_bids_dir, [], [], 0)
% In this case, the function expects appropriate rawdata under path_to_roots_bids_dir
% No other input is necessary, as folders are defined by BIDS standards
%
% (2) Convert DICOM sourcedata into a BIDS compliant dataset and allow selection of files
% that should be used in lead-dbs.
% A typical call would be ea_dataset_import(path_to_input_dir, path_to_output_dir, 1, 1)
%
% input
%   source_dir (mandatory, string or cell containing strings):
%                   (i) a string with the root folder to a BIDS dataset with sourcedata/sub-<label> folder
%                   (ii) a cell with n strings to subject folders, each with a DICOM folder inside containing DICOM raw data
%                        to be converted
%   dest_dir (optional, string or cell containing strings):
%                   (i) string to output folder that is a BIDS dataset rootfolder
%                   (i) a cell with n strings to subject folders, must be the same dimension as source_dir in that case
%   method (integer): DICOM -> Niftii conversion to be used (1 - dcm2niix, 2 - dicm2nii (Matlab), 3 - SPM)
%   dicomimport (boolean integer): if dicomimport is to be triggered or not
% __________________________________________________________________________________
% Copyright (C) 2021 Charite University Medicine Berlin, Movement Disorders Unit
% Johannes Achtzehn

preop_modalities = {'T1w', 'T2w', 'PDw', 'FGATIR'};                         % TODO: get this from prefs
postop_modalities = {'CT', 'ax_MRI', 'sag_MRI', 'cor_MRI'};               % TODO: get this from prefs

if ~iscell(source_dir)   % if single string is passed, convert to cell
    source_dir = {source_dir};
end

if ~iscell(dest_dir) && ~isempty(dest_dir)  % if single string is passed, convert to cell
    dest_dir = {dest_dir};
elseif isempty(dest_dir)     % if dest_dir is emtpy, set it to source_dir
    dest_dir = source_dir;
end

%% import directly from BIDS
if ~dicomimport

    % TODO: also use nifti_to_bids gui?

    if ~contains(source_dir{1}, 'sub-') % root dataset folder has been passed
        rawdata_dir = fullfile(source_dir{1}, 'rawdata');
        lead_derivatives_dir = fullfile(source_dir{1}, 'derivatives', 'leaddbs');

        % find all subjects within the BIDS dataset
        all_files = dir(fullfile(rawdata_dir, 'sub-*'));    % get subjects in dataset root
        dirFlags = [all_files.isdir];                              % get a logical vector that tells which is a directory
        subj_ids = all_files(dirFlags);                         % struct with subject names inside

        fprintf('Found %s subjects and importing directly from BIDS rawdata folder...\n', num2str(numel(subj_ids)))
        for subj_idx = 1:numel(subj_ids)

            fprintf('Importing BIDS data from subject %s...\n', subj_ids(subj_idx).name);

            % preop
            fprintf('\nSearching files for preoperative session...\n');
            anat_files.('preop').anat = find_anat_files_bids(fullfile(rawdata_dir, subj_ids(subj_idx).name, 'ses-preop', 'anat'), preop_modalities);

            % postop
            fprintf('\nSearching files for postoperative session...\n');
            anat_files.('postop').anat = find_anat_files_bids(fullfile(rawdata_dir, subj_ids(subj_idx).name, 'ses-postop', 'anat'), postop_modalities);

            % select which ones to use (if there are multiple)
            anat_files_selected = select_anat_files(anat_files);

            % write into json file
            if ~exist(fullfile(lead_derivatives_dir, subj_ids(subj_idx).name, 'prefs'), 'dir')
                mkdir(fullfile(lead_derivatives_dir, subj_ids(subj_idx).name, 'prefs'));
            end
            savejson('', anat_files_selected, fullfile(lead_derivatives_dir, subj_ids(subj_idx).name, 'prefs', [subj_ids(subj_idx).name, '_desc-rawimages.json']));
        end

    else % only a single rawdata subject folder has been passed
        subjId = regexp(source_dir{1}, ['(?<=rawdata\', filesep, 'sub-).*'], 'match');
        BIDSRoot = regexp(source_dir{1}, ['^.*(?=\', filesep, 'rawdata)'], 'match', 'once');

        fprintf('Importing BIDS data from subject %s...\n', ['sub-', subjId{1}]);

        % preop
        fprintf('\nSearching files for preoperative session...\n');
        anat_files.('preop').anat = find_anat_files_bids(fullfile(source_dir{1}, 'ses-preop', 'anat'), preop_modalities);

        % postop
        fprintf('\nSearching files for postoperative session...\n');
        anat_files.('postop').anat = find_anat_files_bids(fullfile(source_dir{1}, 'ses-postop', 'anat'), postop_modalities);

        % select which ones to use (if there are multiple)
        anat_files_selected = select_anat_files(anat_files);

        % write into json file
        if ~exist(fullfile(BIDSRoot, 'derivatives', 'leaddbs', ['sub-', subjId{1}], 'prefs'), 'dir')
            mkdir(fullfile(BIDSRoot, 'derivatives', 'leaddbs', ['sub-', subjId{1}], 'prefs'));
        end
        savejson('', anat_files_selected, fullfile(BIDSRoot, 'derivatives', 'leaddbs', ['sub-', subjId{1}], 'prefs', [['sub-', subjId{1}], '_desc-rawimages.json']));


    end
    %% import from DICOM
else
    % 1. first check to see if this is already a BIDS root folder or just a
    %     patient folder, create subj_ids for next step

    % case: dataset root is selected and sourcedata is present with
    % according subject folders below it
    if exist(fullfile(source_dir{1}, 'sourcedata'), 'dir') && length(source_dir) == 1
        root_dataset_dir = fullfile(source_dir{1});
        sourcedata_dir = fullfile(source_dir{1}, 'sourcedata');
        lead_derivatives_dir = fullfile(source_dir{1}, 'derivatives', 'leaddbs');
        all_files = dir(fullfile(sourcedata_dir, 'sub-*'));        % get subjects in dataset root
        dirFlags = [all_files.isdir];                              % get a logical vector that tells which is a directory
        subj_ids = all_files(dirFlags);                        % extract only those that are directories (fail-safe)
        fprintf('Found sourcedata for %s subjects and importing DICOMS...\n', num2str(numel(subj_ids)))
        % TODO: check whether there are actually .dcm files present in each
        % subject folder. If not, exclude that subject from the list

        % case: one or more patients selected, but not in BIDS format (e.g. legacy dataset)
    else

        if length(dest_dir) == 1
            fprintf('Found DICOM data for %d subjects, importing DICOMS and creating a BIDS dataset at %s...\n', ...
                length(source_dir), dest_dir{1});
        else
            fprintf('Found DICOM data for %d subjects, importing DICOMS and creating BIDS datasets inside their folders...\n', ...
                length(source_dir));
        end

        subj_ids = struct();
        for subj_idx = 1:length(source_dir)
            if length(dest_dir) == 1                            % this is the case of one destination dir
                dest_dir_subj =  dest_dir{1};
            else
                dest_dir_subj = dest_dir{subj_idx};     % each subj has their own destination dir
            end

            % get subjID from path
            pathparts = strsplit(source_dir{subj_idx}, filesep);
            subjID = pathparts{end};

            % check if sub is already present
            if ~startsWith(subjID, 'sub-')
                subjID = ['sub-', regexprep(subjID, '[\W_]', '')];
            end

            % create folders
            if ~exist(fullfile(dest_dir_subj, 'rawdata', subjID), 'dir')
                mkdir(fullfile(dest_dir_subj, 'rawdata', subjID));
            end
            if ~exist(fullfile(dest_dir_subj, 'sourcedata', subjID), 'dir')
                mkdir(fullfile(dest_dir_subj, 'sourcedata', subjID));
            end

            % copy DICOM data to sourcedata:
            % without DICOM, error when subject only has DICOMS inside
            if ~strcmp(source_dir{subj_idx}, fullfile(dest_dir_subj, 'sourcedata', subjID)) ...
                    && ~exist(fullfile(dest_dir_subj, 'sourcedata', subjID, 'DICOM'), 'dir')
                copyfile(fullfile(source_dir{subj_idx}, 'DICOM'), fullfile(dest_dir_subj, 'sourcedata', subjID));
            end
            subj_ids(subj_idx).name =char(subjID);
        end

        % directory definitions for conversion
        if length(dest_dir) == 1
            lead_derivatives_dir = fullfile(dest_dir{1},'derivatives', 'leaddbs');
            sourcedata_dir = fullfile(dest_dir{1}, 'sourcedata');
            root_dataset_dir = dest_dir{1};
        else
            for subj_idx = 1:length(source_dir)
                lead_derivatives_dir{subj_idx} = fullfile(dest_dir{subj_idx}, 'derivatives', 'leaddbs');
                sourcedata_dir{subj_idx} = fullfile(dest_dir{subj_idx}, 'sourcedata');
                root_dataset_dir{subj_idx} = dest_dir{subj_idx};
            end
        end

    end

    % 2. once subj_ids is created, go through all of the subjects and convert DICOMS
    for subj_idx = 1:numel(subj_ids)

        fprintf('Importing DICOM data from subject %s...\n', subj_ids(subj_idx).name);

        if iscell(sourcedata_dir)   % multiple subjects
            sourcedata_dir_subj = sourcedata_dir{subj_idx};
            root_dataset_dir_subj = root_dataset_dir{subj_idx};
            lead_derivatives_dir_subj = lead_derivatives_dir{subj_idx};
        else    % just one subject
            sourcedata_dir_subj = sourcedata_dir;
            root_dataset_dir_subj = root_dataset_dir;
            lead_derivatives_dir_subj = lead_derivatives_dir;
        end
        dicom_dir = fullfile(sourcedata_dir_subj, subj_ids(subj_idx).name);

        h_wait = waitbar(0, 'Please wait while DICOM files are converted to NIFTI images...');
        niiFiles = ea_dcm_to_nii(dicom_dir, fullfile(dicom_dir, 'tmp'), method); % convert DICOM to NIfTI and get list of files
        close(h_wait);

        % call GUI to select which files should be loaded
        anat_files_selected = ea_nifti_to_bids(niiFiles, root_dataset_dir_subj, subj_ids(subj_idx).name);

        if ~isempty(anat_files_selected)
            % write into json file
            if ~exist(fullfile(lead_derivatives_dir_subj, subj_ids(subj_idx).name, 'prefs'), 'dir')
                mkdir(fullfile(lead_derivatives_dir_subj, subj_ids(subj_idx).name, 'prefs'));
            end
            savejson_struct = struct('filename', fullfile(lead_derivatives_dir_subj, subj_ids(subj_idx).name, 'prefs', [subj_ids(subj_idx).name, '_desc-rawimages.json']), ...
                'singletcell', 0);
            savejson('', anat_files_selected, savejson_struct);
            rmdir(fullfile(dicom_dir, 'tmp'), 's');
        else
            % delete temporary files and folder
            ea_delete(fullfile(dicom_dir, 'tmp'));
            ea_delete(fullfile(root_dataset_dir, 'rawdata', subjID));
        end

        % second option: use lookup table to find files and convert them to BIDS
        % read in lookup table
        %f = fileread(fullfile(ea_getearoot(), 'ext_libs', 'dcm2nii', 'dicom_bids_heuristics.json'));
        %bids_naming_heuristics = jsondecode(f);

        % thirs optpion: heudiconv

        % fourth option: just leave them and user has to manually rename files
        %anat_files.preop = find_anat_files_dicom(tmp_dir, bids_naming_heuristics.preop);

    end

end

    function found_files = find_anat_files_bids(folder, modalities)
        % this function looks in a specific BIDS folder for images that match a specific pattern

        % input
        %   folder -> the folder to look in
        %   modalities -> list of patterns to match files to
        % output
        %   % found_files -> struct for pre-/postop with cells of file names

        found_files = struct();              % store found files in a struct
        for mod_idx = 1:numel(modalities)       % go through the different preop-image types according to prefs and see if we can find some files

            if any(strcmp(modalities{mod_idx}, {'ax_MRI', 'cor_MRI', 'sag_MRI'}))
                % find out which one was found
                mri_dirs = {'ax', 'cor', 'sag'};
                dir_idx = strcmp(modalities{mod_idx}, {'ax_MRI', 'cor_MRI', 'sag_MRI'});
                all_files = dir(fullfile(folder, ['*_acq-', mri_dirs{dir_idx}, '_*', '.nii.gz']));
            else
                all_files = dir(fullfile(folder, ['*_', modalities{mod_idx}, '.nii.gz']));
            end

            switch numel(all_files)
                case 0
                    found_files.(modalities{mod_idx}) = {};
                otherwise
                    found_files_list = {};
                    for name_idx = 1:length(all_files)
                        fprintf('Found file %s for modality %s\n', all_files(name_idx).name, modalities{mod_idx});

                        % remove extension from filename
                        dot_idx = strfind(all_files(name_idx).name, '.');

                        if length(dot_idx) > 1  % for .nii.gz, we have 2 dots...
                            dot_idx = dot_idx(1);
                        end

                        fname = all_files(name_idx).name(1:dot_idx - 1);
                        found_files_list{name_idx, 1} = fname;
                    end
                    found_files.(modalities{mod_idx}) = found_files_list;
            end
        end
    end

    function anat_files = select_anat_files(anat_files)
        % this function checks if multiple files have been found for each pre-/postop modality and prompts the user to select on of them

        % input
        %   anat_files -> struct for pre-/postop image file names and cells for each filename
        % output
        %   anat_files -> same as input, but with fields removed that have no filenames. Multiple entries are also deleted

        % preop
        for mod_idx = 1:length(fieldnames(anat_files.preop.anat))

            if numel(anat_files.preop.anat.(preop_modalities{mod_idx})) > 1         % if more than one is found
                file_list = anat_files.preop.anat.(preop_modalities{mod_idx});

                % prompt the user to select a preop image, here we only allow one selection
                [index, tf] = listdlg('PromptString', {sprintf('Multiple preop files for modality %s found',preop_modalities{mod_idx}), ...
                    'Please select on file from the list'}, 'Name', 'Select preop images', 'SelectionMode', 'single', 'ListSize', [300, 500], 'ListString', file_list);
                anat_files.preop.anat.(preop_modalities{mod_idx}) = file_list{index, 1};

            elseif numel(anat_files.preop.anat.(preop_modalities{mod_idx})) == 0   % if none is found
                anat_files.preop.anat = rmfield(anat_files.preop.anat, preop_modalities{mod_idx});

            else    % if only one convert (for better json readability)
                anat_files.preop.anat.(preop_modalities{mod_idx}) = cell2mat(anat_files.preop.anat.(preop_modalities{mod_idx}));
            end
        end

        % postop
        for mod_idx = 1:length(fieldnames(anat_files.postop.anat))

            if numel(anat_files.postop.anat.(postop_modalities{mod_idx})) > 1   % if more than one is found
                file_list = anat_files.postop.anat.(postop_modalities{mod_idx});

                % prompt the user to select postop images, for everything other than CT, allow multiple selections (MRI)
                [index, tf] = listdlg('PromptString', {sprintf('Multiple postop files for modality %s found',postop_modalities{mod_idx}), ...
                    'Please select on file from the list'}, 'Name', 'Select postop images', 'SelectionMode', 'single', ...
                    'ListSize', [300, 500], 'ListString', file_list);

                anat_files.postop.anat.(postop_modalities{mod_idx}) = file_list{index, 1};

            elseif numel(anat_files.postop.anat.(postop_modalities{mod_idx})) == 0   % if none is found
                anat_files.postop.anat = rmfield(anat_files.postop.anat, postop_modalities{mod_idx});

            else     % if only one convert (for better json readability)
                anat_files.postop.anat.(postop_modalities{mod_idx}) = cell2mat(anat_files.postop.anat.(postop_modalities{mod_idx}));
            end
        end
    end

end
