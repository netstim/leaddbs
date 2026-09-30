function options = ea_prepare_dti_bids(options)
% Copy DWI files from BIDS rawdata to derivatives with BIDS-compliant names

directory = [options.root, options.patientname, filesep];

% Create preprocessing/dwi directory if needed
derivDwiDir = fullfile(directory, 'preprocessing', 'dwi');
if ~isfolder(derivDwiDir)
    mkdir(derivDwiDir);
end

% Find rawdata directory
try
    % Extract dataset root from derivatives path
    datasetRoot = regexp(directory, ['^.*(?=\', filesep, 'derivatives)'], 'match', 'once');
    subjId = regexp(options.patientname, '(?<=sub-).+', 'match', 'once');
    if isempty(subjId)
        subjId = options.patientname;
    end
    
    rawDataDir = fullfile(datasetRoot, 'rawdata', ['sub-', subjId]);
    
    % Find DWI files in rawdata (search recursively)
    % Accept both BIDS standard '_dwi' and legacy '_DTI' suffixes
    rawDwiFiles = [dir(fullfile(rawDataDir, '**', '*_dwi.nii.gz')); ...
                   dir(fullfile(rawDataDir, '**', '*_DTI.nii.gz'))];
    if isempty(rawDwiFiles)
        rawDwiFiles = [dir(fullfile(rawDataDir, '**', '*_dwi.nii')); ...
                       dir(fullfile(rawDataDir, '**', '*_DTI.nii'))];
    end

    if ~isempty(rawDwiFiles)
        % Use first DWI run found
        rawDwiPath = fullfile(rawDwiFiles(1).folder, rawDwiFiles(1).name);
        [~, dwiBaseName, dwiExt] = fileparts(rawDwiPath);
        if strcmp(dwiExt, '.gz')
            [~, dwiBaseName] = fileparts(dwiBaseName);
        end

        % Normalize suffix to '_dwi' regardless of source naming (_DTI -> _dwi)
        dwiBaseNameNorm = regexprep(dwiBaseName, '_DTI$', '_dwi');

        % Target BIDS-compliant file in derivatives/preprocessing/dwi/
        targetDwi = fullfile(derivDwiDir, [dwiBaseNameNorm, '.nii']);

        if ~exist(targetDwi, 'file')
            if strcmp(dwiExt, '.gz')
                % Gunzip to temp name then rename if normalization was needed
                gunzip(rawDwiPath, derivDwiDir);
                if ~strcmp(dwiBaseName, dwiBaseNameNorm)
                    movefile(fullfile(derivDwiDir, [dwiBaseName, '.nii']), targetDwi);
                end
            else
                % Copy with normalized name
                copyfile(rawDwiPath, targetDwi);
            end
        end

        % Copy .bval and .bvec (source uses original name, target uses normalized)
        rawDwiDir = rawDwiFiles(1).folder;
        rawBval = fullfile(rawDwiDir, [dwiBaseName, '.bval']);
        rawBvec = fullfile(rawDwiDir, [dwiBaseName, '.bvec']);
        targetBval = fullfile(derivDwiDir, [dwiBaseNameNorm, '.bval']);
        targetBvec = fullfile(derivDwiDir, [dwiBaseNameNorm, '.bvec']);

        if exist(rawBval, 'file') && ~exist(targetBval, 'file')
            copyfile(rawBval, targetBval);
        end
        if exist(rawBvec, 'file') && ~exist(targetBvec, 'file')
            copyfile(rawBvec, targetBvec);
        end

        % Update options.prefs to point to BIDS paths in preprocessing/dwi
        options.prefs.dti = fullfile('preprocessing', 'dwi', [dwiBaseNameNorm, '.nii']);
        options.prefs.bval = fullfile('preprocessing', 'dwi', [dwiBaseNameNorm, '.bval']);
        options.prefs.bvec = fullfile('preprocessing', 'dwi', [dwiBaseNameNorm, '.bvec']);
        options.prefs.b0 = fullfile('preprocessing', 'dwi', ['sub-', subjId, '_ses-preop_b0.nii']);
        options.prefs.fa = fullfile('preprocessing', 'dwi', ['sub-', subjId, '_ses-preop_fa.nii']);
        
        disp(['DWI files prepared: ', targetDwi]);
    else
        warning('No DWI files found in rawdata for subject %s', subjId);
    end
catch ME
    warning('Failed to prepare DWI files: %s', ME.message);
end

