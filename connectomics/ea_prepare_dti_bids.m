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
    rawDwiFiles = dir(fullfile(rawDataDir, '**', '*_dwi.nii.gz'));
    if isempty(rawDwiFiles)
        rawDwiFiles = dir(fullfile(rawDataDir, '**', '*_dwi.nii'));
    end
    
    if ~isempty(rawDwiFiles)
        % Use first DWI run found
        rawDwiPath = fullfile(rawDwiFiles(1).folder, rawDwiFiles(1).name);
        [~, dwiBaseName, dwiExt] = fileparts(rawDwiPath);
        if strcmp(dwiExt, '.gz')
            [~, dwiBaseName] = fileparts(dwiBaseName);
        end
        
        % Target BIDS-compliant file in derivatives/preprocessing/dwi/
        targetDwi = fullfile(derivDwiDir, [dwiBaseName, '.nii']);
        
        if ~exist(targetDwi, 'file')
            disp(['Copying DWI from rawdata: ', dwiBaseName]);
            if strcmp(dwiExt, '.gz')
                % Gunzip
                gunzip(rawDwiPath, derivDwiDir);
            else
                % Copy
                copyfile(rawDwiPath, targetDwi);
            end
        end
        
        % Copy .bval and .bvec
        rawDwiDir = rawDwiFiles(1).folder;
        rawDwiName = dwiBaseName;  % Already cleaned above
        rawBval = fullfile(rawDwiDir, [rawDwiName, '.bval']);
        rawBvec = fullfile(rawDwiDir, [rawDwiName, '.bvec']);
        targetBval = fullfile(derivDwiDir, [dwiBaseName, '.bval']);
        targetBvec = fullfile(derivDwiDir, [dwiBaseName, '.bvec']);
        
        if exist(rawBval, 'file') && ~exist(targetBval, 'file')
            copyfile(rawBval, targetBval);
        end
        if exist(rawBvec, 'file') && ~exist(targetBvec, 'file')
            copyfile(rawBvec, targetBvec);
        end
        
        % Update options.prefs to point to BIDS paths in preprocessing/dwi
        options.prefs.dti = fullfile('preprocessing', 'dwi', [dwiBaseName, '.nii']);
        options.prefs.bval = fullfile('preprocessing', 'dwi', [dwiBaseName, '.bval']);
        options.prefs.bvec = fullfile('preprocessing', 'dwi', [dwiBaseName, '.bvec']);
        options.prefs.b0 = fullfile('preprocessing', 'dwi', [dwiBaseName, '_b0.nii']);
        options.prefs.fa = fullfile('preprocessing', 'dwi', [dwiBaseName, '_fa.nii']);
        
        disp(['DWI files prepared: ', targetDwi]);
    else
        warning('No DWI files found in rawdata for subject %s', subjId);
    end
catch ME
    warning('Failed to prepare DWI files: %s', ME.message);
end

