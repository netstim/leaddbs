function ea_dicom_import(options)
% This function converts DICOM files from the sourcedata folder of your dataset into
% a BIDS-conform rawdata folder and specified which files are to be used by lead-dbs
% __________________________________________________________________________________
% Copyright (C) 20121 Charite University Medicine Berlin, Movement Disorders Unit
% Johannes Achtzehn & Andreas Horn

preop_modalities = {'T1w', 'T2w', 'PDw', 'FGATIR'};                 % TODO: get this from prefs
postop_modalities = {'CT', 'ax_MR', 'sag_MR', 'cor_MR'};            % TODO: get this from prefs

% TODO: 1. GUI to select which files, 2. first get how many sessions, then iterate over those (currently only pre- and postop)
% TODO: 3. handle no files gracefully

%% import directly from BIDS
if isfield(options.dicomimp,'method') && options.dicomimp.method == 4
    
    rawdata_dir = fullfile(options.root, options.patientname, 'rawdata');   % TODO: dataset fetcher
    lead_derivatives_dir = fullfile(options.root, options.patientname, 'derivatives', 'leaddbs');
    
    all_files = dir(fullfile(rawdata_dir, 'sub-*'));    % get subjects in dataset root TODO: dataset fetcher
    dirFlags = [all_files.isdir];                       % get a logical vector that tells which is a directory
    subj_folders = all_files(dirFlags);                 % extract only those that are directories (fail-safe)
    
    fprintf('Found %s subjects and importing directly from BIDS rawdata folder...\n', num2str(numel(subj_folders)))
    for subj_idx = 1:numel(subj_folders)
        
        fprintf('Importing BIDS data from subject %s...\n', subj_folders(subj_idx).name);
        
        % preop
        fprintf('\nSearching files for preoperative session...\n');
        anat_files.('preop') = find_anat_files_bids(fullfile(rawdata_dir, subj_folders(subj_idx).name, 'ses-preop', 'anat'), preop_modalities);
        
        % postop
        fprintf('\nSearching files for postoperative session...\n');
        anat_files.('postop') = find_anat_files_bids(fullfile(rawdata_dir, subj_folders(subj_idx).name, 'ses-postop', 'anat'), postop_modalities);
        
        % select which ones to use (if there are multiple)
        anat_files_selected = select_anat_files(anat_files);
        
        % write into json file
        if ~exist(fullfile(lead_derivatives_dir, subj_folders(subj_idx).name, 'prefs'), 'dir')
            mkdir(fullfile(lead_derivatives_dir, subj_folders(subj_idx).name, 'prefs'));
        end
        savejson('', anat_files_selected, fullfile(lead_derivatives_dir, subj_folders(subj_idx).name, 'prefs', [subj_folders(subj_idx).name, '_desc-rawimages.json']));
    end
    
    %% import from DICOM
else
    sourcedata_dir = fullfile(options.root, options.patientname, 'sourcedata');     % TODO: dataset fetcher
    
    all_files = dir(fullfile(sourcedata_dir, 'sub-*'));  % get subjects in dataset root TODO: dataset fetcher
    dirFlags = [all_files.isdir];                        % get a logical vector that tells which is a directory
    subj_folders = all_files(dirFlags);                  % extract only those that are directories (fail-safe)
    
    for subj_idx = 1:numel(subj_folders)
        
        fprintf('Importing DICOM data from subject %s...\n', subj_folders(subj_idx).name);
        
        dicom_dir = fullfile(options.root, options.patientname, 'sourcedata', subj_folders(subj_idx).name);
        tmp_dir = fullfile(options.root, options.patientname, 'sourcedata', subj_folders(subj_idx).name, 'tmp'); % tmp dir for temporary storage of niftis
        
%         switch options.dicomimp.method
%             case 1 % dcm2niix
%                 ea_dcm2niix(dicom_dir, tmp_dir);
%                 ea_methods(options, ['DICOM images were converted to the '...
%                     'NIfTI file format using dcm2niix (see https://github.com/rordenlab/dcm2niix).']);
%             case 2 % dicm2nii
%                 ea_dicm2nii(dicom_dir, tmp_dir);
%                 ea_methods(options, ['DICOM images were converted to the '...
%                     'NIfTI file format using dicm2nii (see https://github.com/xiangruili/dicm2nii).']);
%             case 3 % SPM
%                 ea_spm_dicom_import(dicom_dir, tmp_dir);
%                 ea_methods(options, ['DICOM images were converted to the '...
%                     'NIfTI file format using SPM.']);
%         end
%         
        % first option: rename files with the help of a GUI
        [~, niiFiles] = fileparts(ea_regexpdir(tmp_dir, '\.nii.gz$', 0));
        ea_select_nii_import(niiFiles, tmp_dir, dicom_dir);
        disp('done');
        % second option: use lookup table to find files and convert them to BIDS
        % read in lookup table
        %f = fileread(fullfile(ea_getearoot(), 'ext_libs', 'dcm2nii', 'dicom_bids_heuristics.json'));
        %bids_naming_heuristics = jsondecode(f);
        
        % thirs optpion: heudiconv
        
        % fourth option: just leave them and user has to manually rename files
        %anat_files.preop = find_anat_files_dicom(tmp_dir, bids_naming_heuristics.preop);
        
    end
    
    % methods call
    switch options.dicomimp.method
        case 1 % dcm2niix
            ea_methods(options, ['DICOM images were converted to the '...
                'NIfTI file format using dcm2niix (see https://github.com/rordenlab/dcm2niix).']);
        case 2 % dicm2nii
            ea_methods(options, ['DICOM images were converted to the '...
                'NIfTI file format using dicm2nii (see https://github.com/xiangruili/dicm2nii).']);
        case 3 % SPM
            ea_methods(options, ['DICOM images were converted to the '...
                'NIfTI file format using SPM.']);
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

    function found_files = find_anat_files_dicom(folder, patterns)
        found_files = struct();
        
        for mod_idx = 1:numel(fields(patterns))
            patterns_mod = patterns.(patterns)
            
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
