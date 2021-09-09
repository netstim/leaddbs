function ea_dataset_import(source_dir, dest_dir, method, dicomimport)
% This function converts datasets from the sourcedata folder of your dataset into
% a BIDS-conform rawdata folder and specified which files are to be used by lead-dbs
% __________________________________________________________________________________
% Copyright (C) 20121 Charite University Medicine Berlin, Movement Disorders Unit
% Johannes Achtzehn & Andreas Horn

preop_modalities = {'T1w', 'T2w', 'PDw', 'FGATIR'};                 % TODO: get this from prefs
postop_modalities = {'CT', 'ax_MR', 'sag_MR', 'cor_MR'};            % TODO: get this from prefs

if ~iscell(source_dir)
    source_dir = {source_dir};
end

if ~iscell(dest_dir)
    dest_dir = {dest_dir};
elseif isempty(dest_dir)
    dest_dir = source_dir;
end

%% import directly from BIDS
if ~dicomimport
    
    % TODO: also use dicom_to_bids gui?
    
    % if source_dir is passed as cell, convert it
    
    rawdata_dir = fullfile(source_dir{1}, 'rawdata');
    lead_derivatives_dir = fullfile(source_dir{1}, 'derivatives', 'leaddbs');
    
    % before running, lets
    all_files = dir(fullfile(rawdata_dir, 'sub-*'));    % get subjects in dataset root
    dirFlags = [all_files.isdir];                       % get a logical vector that tells which is a directory
    subj_ids = all_files(dirFlags);                 % extract only those that are directories (fail-safe)
    
    fprintf('Found %s subjects and importing directly from BIDS rawdata folder...\n', num2str(numel(subj_ids)))
    for subj_idx = 1:numel(subj_ids)
        
        fprintf('Importing BIDS data from subject %s...\n', subj_ids(subj_idx).name);
        
        % preop
        fprintf('\nSearching files for preoperative session...\n');
        anat_files.('preop') = find_anat_files_bids(fullfile(rawdata_dir, subj_ids(subj_idx).name, 'ses-preop', 'anat'), preop_modalities);
        
        % postop
        fprintf('\nSearching files for postoperative session...\n');
        anat_files.('postop') = find_anat_files_bids(fullfile(rawdata_dir, subj_ids(subj_idx).name, 'ses-postop', 'anat'), postop_modalities);
        
        % select which ones to use (if there are multiple)
        anat_files_selected = select_anat_files(anat_files);
        
        % write into json file
        if ~exist(fullfile(lead_derivatives_dir, subj_ids(subj_idx).name, 'prefs'), 'dir')
            mkdir(fullfile(lead_derivatives_dir, subj_ids(subj_idx).name, 'prefs'));
        end
        savejson('', anat_files_selected, fullfile(lead_derivatives_dir, subj_ids(subj_idx).name, 'prefs', [subj_ids(subj_idx).name, '_desc-rawimages.json']));
    end
    
    %% import from DICOM
else
    % check to see if this is already a BIDS root folder or just a
    % patient folder, create subj_ids for next step
    
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
        
    % case: one or more patients selected, but not in BIDS format
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
            
            pathparts = strsplit(source_dir{subj_idx}, filesep);
            subjID = char(pathparts(end));
            
            % check if sub is already present
            if ~strcmp(subjID(1:4), 'sub-')
                subjID = ['sub-', subjID];
            end

            % create folders
            if ~exist(fullfile(dest_dir_subj, 'rawdata', subjID), 'dir')
                mkdir(fullfile(dest_dir_subj, 'rawdata', subjID));
            end
            if ~exist(fullfile(dest_dir_subj, 'sourcedata', subjID), 'dir')
                mkdir(fullfile(dest_dir_subj, 'sourcedata', subjID));
            end
            
            copyfile(fullfile(source_dir{subj_idx}, 'DICOM'), fullfile(dest_dir_subj, 'sourcedata', subjID));
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
    
    % now go through all of the subjects and convert DICOMS
    for subj_idx = 1:numel(subj_ids)
        
        fprintf('Importing DICOM data from subject %s...\n', subj_ids(subj_idx).name);
        
        if iscell(sourcedata_dir)
            sourcedata_dir_subj = sourcedata_dir{subj_idx};
            root_dataset_dir_subj = root_dataset_dir{subj_idx};
            lead_derivatives_dir_subj = lead_derivatives_dir{subj_idx};
        else
            sourcedata_dir_subj = sourcedata_dir;
            root_dataset_dir_subj = root_dataset_dir;
            lead_derivatives_dir_subj = lead_derivatives_dir;
        end
        dicom_dir = fullfile(sourcedata_dir_subj, subj_ids(subj_idx).name);
        
        % convert DICOMS to .nii files and get list of files
        niiFiles = ea_dcm_to_nii(method, dicom_dir);
        
        % call GUI to select which files should be loaded
        anat_files_selected = ea_dicom_to_bids(subj_ids(subj_idx).name, niiFiles, root_dataset_dir_subj);
        
        if ~isempty(anat_files_selected)
            % write into json file
            if ~exist(fullfile(lead_derivatives_dir_subj, subj_ids(subj_idx).name, 'prefs'), 'dir')
                mkdir(fullfile(lead_derivatives_dir_subj, subj_ids(subj_idx).name, 'prefs'));
            end
            savejson('', anat_files_selected, fullfile(lead_derivatives_dir_subj, subj_ids(subj_idx).name, 'prefs', [subj_ids(subj_idx).name, '_desc-rawimages.json']));
        end
        
        % second option: use lookup table to find files and convert them to BIDS
        % read in lookup table
        %f = fileread(fullfile(ea_getearoot(), 'ext_libs', 'dcm2nii', 'dicom_bids_heuristics.json'));
        %bids_naming_heuristics = jsondecode(f);
        
        % thirs optpion: heudiconv
        
        % fourth option: just leave them and user has to manually rename files
        %anat_files.preop = find_anat_files_dicom(tmp_dir, bids_naming_heuristics.preop);
        
        rmdir(fullfile(dicom_dir, 'tmp'), 's');
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
            
            if any(strcmp(modalities{mod_idx}, {'ax_MR', 'cor_MR', 'sag_MR'}))
                % find out which one was found
                mri_dirs = {'ax', 'cor', 'sag'};
                dir_idx = strcmp(modalities{mod_idx}, {'ax_MR', 'cor_MR', 'sag_MR'});
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
                        found_files_list{name_idx, 1} = all_files(name_idx).name;
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
        for mod_idx = 1:length(fieldnames(anat_files.preop))
            
            if numel(anat_files.preop.(preop_modalities{mod_idx})) > 1         % if more than one is found
                file_list = anat_files.preop.(preop_modalities{mod_idx});
                
                % prompt the user to select a preop image, here we only allow one selection
                [index, tf] = listdlg('PromptString', {sprintf('Multiple preop files for modality %s found',preop_modalities{mod_idx}), ...
                    'Please select on file from the list'}, 'Name', 'Select preop images', 'SelectionMode', 'single', 'ListSize', [300, 500], 'ListString', file_list);
                anat_files.preop.(preop_modalities{mod_idx}) = file_list{index, 1};
                
            elseif numel(anat_files.preop.(preop_modalities{mod_idx})) == 0   % if none is found
                anat_files.preop = rmfield(anat_files.preop, preop_modalities{mod_idx});
                
            else    % if only one convert (for better json readability)
                anat_files.preop.(preop_modalities{mod_idx}) = cell2mat(anat_files.preop.(preop_modalities{mod_idx}));
            end
        end
        
        % postop
        for mod_idx = 1:length(fieldnames(anat_files.postop))
            
            if numel(anat_files.postop.(postop_modalities{mod_idx})) > 1   % if more than one is found
                file_list = anat_files.postop.(postop_modalities{mod_idx});
                
                % prompt the user to select postop images, for everything other than CT, allow multiple selections (MRI)
                [index, tf] = listdlg('PromptString', {sprintf('Multiple postop files for modality %s found',postop_modalities{mod_idx}), ...
                    'Please select on file from the list'}, 'Name', 'Select postop images', 'SelectionMode', 'single', ...
                    'ListSize', [300, 500], 'ListString', file_list);
                
                anat_files.postop.(postop_modalities{mod_idx}) = file_list{index, 1};
                
            elseif numel(anat_files.postop.(postop_modalities{mod_idx})) == 0   % if none is found
                anat_files.postop = rmfield(anat_files.postop, postop_modalities{mod_idx});
                
            else     % if only one convert (for better json readability)
                anat_files.postop.(postop_modalities{mod_idx}) = cell2mat(anat_files.postop.(postop_modalities{mod_idx}));
            end
        end
    end

end
