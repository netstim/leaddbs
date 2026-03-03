function options = ea_getptopts(directory,options)
% Generate minimal options struct from a patient directory.

if isempty(directory)
    directory = pwd;
end

directory = GetFullPath(directory);

% BIDS FIX: Remove trailing filesep for consistent path handling
if endsWith(directory, filesep)
    directory = directory(1:end-1);
end

options.earoot = ea_getearoot;
options.native = 0;

if contains(directory, ['derivatives', filesep, 'leaddbs'])
    % Construct BIDS class, get subjId
    BIDSRoot = regexp(directory, ['^.*(?=\', filesep, 'derivatives)'], 'match', 'once');
    bids = BIDSFetcher(BIDSRoot);
    subjId = regexp(directory, ['(?<=leaddbs\', filesep, 'sub-)[^\', filesep, ']+'], 'match', 'once');

    % Set prefs (will be overridden with BIDS paths later)
    options.prefs = ea_prefs(['sub-', subjId]);

    % Set modality field
    uiPrefsFile = bids.getPrefs(subjId, 'uiprefs', 'mat');
    if isfile(uiPrefsFile)
        uiprefs = load(uiPrefsFile);
        options.modality = uiprefs.modality;
        options.sides = uiprefs.sides;
    else
        options.modality = bids.settings.preferMRCT;
        options.sides = [1,2];
    end

    % Set subj BIDS struct
    try
        options.subj = bids.getSubj(subjId, options.modality);

        % --- Register submarine folder (if it exists) ---
        submarineBase = fullfile(options.subj.subjDir, 'submarine');

        if exist(submarineBase, 'dir')

            % Ensure submarine struct exists
            if ~isfield(options.subj, 'submarine') || isempty(options.subj.submarine)
                options.subj.submarine = struct();
            end

            % Store base submarine directory
            options.subj.submarineDir = submarineBase;

            % Ensure atlases container exists
            if ~isfield(options.subj.submarine, 'atlases') || isempty(options.subj.submarine.atlases)
                options.subj.submarine.atlases = struct();
            end

            % List atlas subfolders
            d = dir(submarineBase);
            d = d([d.isdir] & ~ismember({d.name}, {'.','..'}));

            for k = 1:numel(d)
                atlas_name  = d(k).name;
                atlasField = matlab.lang.makeValidName(atlas_name);

                % Register atlas directory
                options.subj.submarine.atlases.(atlasField).dir = ...
                    fullfile(submarineBase, atlas_name);
            end
        end


        % Set primary template
        subjAnchor = regexprep(options.subj.AnchorModality, '[^\W_]+_', '');
        if ismember(subjAnchor, fieldnames(bids.spacedef.norm_mapping))
            options.primarytemplate = bids.spacedef.norm_mapping.(subjAnchor);
        else
            options.primarytemplate = bids.spacedef.misfit_template;
        end

        % Set elmodel and elspec
        recon = bids.getRecon(subjId);
    catch
        options.subj.subjId = subjId;
        options.subj.subjDir = directory;

        reconFile = ea_regexpdir(fullfile(directory, 'reconstruction'), '_desc-reconstruction\.mat$', 0);
        if ~isempty(reconFile)
            recon.recon = reconFile{1};
            options.subj.recon = recon;
        else
            recon.recon = '';
        end

        options.subj.stimDir = fullfile(directory, 'stimulations');

        statsFile = ea_regexpdir(directory, '_desc-stats\.mat$', 0);
        if ~isempty(statsFile)
            options.subj.stats = statsFile{1};
        end
    end

    % Set elmodel and elspec
    if isfile(recon.recon)
        load(recon.recon);
        try
            options.elmodel = ea_get_first_notempty_elmodel(reco.props);
        catch
            options.elmodel = 'Medtronic 3389';
        end
        options = ea_resolve_elspec(options);
    end

    options.bids=bids; % store bidsfetcher obj within options.

    % Provide classic fields for Lead-Connectome compatibility
    try
        options.root = [fileparts(options.subj.subjDir), filesep];
        [~, options.patientname] = fileparts(options.subj.subjDir);
    catch
    end

    % Point to preprocessed BIDS files if they exist
    try
        preprocFuncDir = fullfile(options.subj.subjDir, 'preprocessing', 'func');
        preprocDwiDir = fullfile(options.subj.subjDir, 'preprocessing', 'dwi');
        preprocAnatDir = fullfile(options.subj.subjDir, 'preprocessing', 'anat');

        % --- Register preprocessing dirs in options.subj.preproc ---
        if ~isfield(options.subj,'preproc') || isempty(options.subj.preproc)
            options.subj.preproc = struct();
        end

        % keep any existing anat field (often set by bids.getSubj)
        if ~isfield(options.subj.preproc,'anat')
            options.subj.preproc.anat = struct();
        end
        options.subj.preproc.anat.dir = preprocAnatDir;

        % add missing ones
        options.subj.preproc.dwi  = struct('dir', preprocDwiDir);
        options.subj.preproc.func = struct('dir', preprocFuncDir);


        % rs-fMRI: Copy from rawdata to preprocessing/func if not present
        if ~isfolder(preprocFuncDir)
            mkdir(preprocFuncDir);
        end

        % Search for BOLD in preprocessing/func first
        boldFiles = dir(fullfile(preprocFuncDir, '*_bold.nii'));
        if isempty(boldFiles)
            boldFiles = dir(fullfile(preprocFuncDir, '*_bold.nii.gz'));
        end

        % If not found, copy from rawdata
        if isempty(boldFiles)
            % Navigate from derivatives/leaddbs/sub-XXX to rawdata/sub-XXX
            derivLeaddbsDir = fileparts(options.subj.subjDir); % .../derivatives/leaddbs
            derivDir = fileparts(derivLeaddbsDir); % .../derivatives
            datasetRoot = fileparts(derivDir); % dataset root
            rawDataDir = fullfile(datasetRoot, 'rawdata', ['sub-', options.subj.subjId]);
            rawBoldFiles = dir(fullfile(rawDataDir, '**', '*_bold.nii.gz'));
            if isempty(rawBoldFiles)
                rawBoldFiles = dir(fullfile(rawDataDir, '**', '*_bold.nii'));
            end

            if ~isempty(rawBoldFiles)
                % Copy first BOLD file
                srcBold = fullfile(rawBoldFiles(1).folder, rawBoldFiles(1).name);
                [~, boldBaseName, boldExt] = fileparts(rawBoldFiles(1).name);
                if strcmp(boldExt, '.gz')
                    [~, boldBaseName] = fileparts(boldBaseName); % remove .nii from .nii.gz
                end
                targetBold = fullfile(preprocFuncDir, [boldBaseName, '.nii']);

                if ~exist(targetBold, 'file')
                    fprintf('Copying rs-fMRI data from rawdata to preprocessing/func...\n');
                    if endsWith(srcBold, '.gz')
                        % Gunzip and copy
                        copyfile(srcBold, [targetBold, '.gz']);
                        gunzip([targetBold, '.gz']);
                        delete([targetBold, '.gz']);
                    else
                        copyfile(srcBold, targetBold);
                    end
                    fprintf('rs-fMRI data copied successfully.\n');
                end

                boldFiles = dir(fullfile(preprocFuncDir, '*_bold.nii'));
            end
        end

        % Set rest and pprest paths
        if ~isempty(boldFiles)
            fullBoldPath = fullfile(boldFiles(1).folder, boldFiles(1).name);
            options.prefs.rest = strrep(fullBoldPath, [options.subj.subjDir, filesep], '');
            % pprest is the smoothed+realigned version (sr prefix)
            [fDir, fName, fExt] = fileparts(fullfile(boldFiles(1).folder, boldFiles(1).name));
            srFile = fullfile(fDir, ['sr', fName, fExt]);
            if exist(srFile, 'file')
                options.prefs.pprest = strrep(srFile, [options.subj.subjDir, filesep], '');
            end
        end

        % DWI: Copy from rawdata to preprocessing/dwi with BIDS-compliant names if not present
        derivDwiDir = fullfile(options.subj.subjDir, 'preprocessing', 'dwi');
        if ~isfolder(derivDwiDir)
            mkdir(derivDwiDir);
        end

        % First check if DWI already exists in preprocessing/dwi
        dwiFiles = dir(fullfile(derivDwiDir, '*_dwi.nii'));

        if isempty(dwiFiles)
            % If not found, look in rawdata (search recursively)
            rawDataDir = fullfile(fileparts(fileparts(options.subj.subjDir)), 'rawdata', ['sub-', options.subj.subjId]);
            rawDwiFiles = dir(fullfile(rawDataDir, '**', '*_dwi.nii.gz'));
            if isempty(rawDwiFiles)
                rawDwiFiles = dir(fullfile(rawDataDir, '**', '*_dwi.nii'));
            end
        else
            % DWI already in preprocessing - use it
            rawDwiFiles = [];
        end

        if ~isempty(dwiFiles)
            % DWI already exists in preprocessing/dwi - set paths directly
            dwiPath = fullfile(dwiFiles(1).folder, dwiFiles(1).name);
            [~, dwiBaseName] = fileparts(dwiPath);

            options.prefs.dti = strrep(dwiPath, [options.subj.subjDir, filesep], '');
            options.prefs.bval = fullfile('preprocessing', 'dwi', [dwiBaseName, '.bval']);
            options.prefs.bvec = fullfile('preprocessing', 'dwi', [dwiBaseName, '.bvec']);
            options.prefs.b0 = fullfile('preprocessing', 'dwi', [dwiBaseName, '_b0.nii']);
            options.prefs.fa = fullfile('preprocessing', 'dwi', [dwiBaseName, '_fa.nii']);
            options.prefs.FTR_unnormalized = fullfile('connectomics', 'dMRI', 'FTR.mat');
            
            options.prefs.FTR_normalized   = fullfile('connectomics', 'dMRI', 'FTR_normalized.mat');

            % --- Store in options.subj.preproc.dwi as well ---
            options.subj.preproc.dwi.dwi = dwiPath;

            % If you want absolute bval/bvec too
            options.subj.preproc.dwi.bval = fullfile(options.subj.subjDir, options.prefs.bval);
            options.subj.preproc.dwi.bvec = fullfile(options.subj.subjDir, options.prefs.bvec);

            % Derived files (may not exist yet, but good to define expected locations)
            options.subj.preproc.dwi.b0 = fullfile(options.subj.subjDir, options.prefs.b0);
            options.subj.preproc.dwi.fa = fullfile(options.subj.subjDir, options.prefs.fa);

        elseif ~isempty(rawDwiFiles)
            % Use first DWI run found
            rawDwiPath = fullfile(rawDwiFiles(1).folder, rawDwiFiles(1).name);
            [~, dwiBaseName, dwiExt] = fileparts(rawDwiPath);
            if strcmp(dwiExt, '.gz')
                [~, dwiBaseName] = fileparts(dwiBaseName);
            end

            % Target BIDS-compliant file in derivatives/preprocessing/dwi/
            targetDwi = fullfile(derivDwiDir, [dwiBaseName, '.nii']);

            if ~exist(targetDwi, 'file')
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

            % Set prefs to BIDS paths (relative to subject dir)
            if exist(targetDwi, 'file')
                options.prefs.dti = strrep(targetDwi, [options.subj.subjDir, filesep], '');
                options.prefs.bval = strrep(targetBval, [options.subj.subjDir, filesep], '');
                options.prefs.bvec = strrep(targetBvec, [options.subj.subjDir, filesep], '');
                % b0 and fa will be generated in same directory
                options.prefs.b0 = fullfile('preprocessing', 'dwi', [dwiBaseName, '_b0.nii']);
                options.prefs.fa = fullfile('preprocessing', 'dwi', [dwiBaseName, '_fa.nii']);
                % Fiber tracking output (BIDS: stored in connectomics/dMRI/)
                options.prefs.FTR_unnormalized = fullfile('connectomics', 'dMRI', 'FTR.mat');
                options.prefs.FTR_normalized   = fullfile('connectomics', 'dMRI', 'FTR_normalized.mat');
            end
        end

        % Anatomical anchor: Copy from rawdata to preprocessing/anat if not present
        % Ensure preprocessing/anat exists
        if ~isfolder(preprocAnatDir)
            mkdir(preprocAnatDir);
        end

        % Search for anatomical files in preprocessing/anat first
        % EXCLUDE SPM segmentation files (c1, c2, c3 prefix) and preprocessed files
        anatFiles = dir(fullfile(preprocAnatDir, '*_T1w.nii'));
        if isempty(anatFiles)
            anatFiles = dir(fullfile(preprocAnatDir, '*_T2w.nii'));
        end

        % Filter out preprocessed versions immediately
        if ~isempty(anatFiles)
            validAnatFiles = [];
            for i = 1:length(anatFiles)
                fname = anatFiles(i).name;
                if ~startsWith(fname, {'c', 'r', 'sr', 'mean'})
                    validAnatFiles = [validAnatFiles; anatFiles(i)];
                end
            end
            anatFiles = validAnatFiles;
        end

        % If not found, copy from rawdata
        if isempty(anatFiles)
            fprintf('No anatomical files found in preprocessing/anat, checking rawdata...\n');
            % Navigate from derivatives/leaddbs/sub-XXX to rawdata/sub-XXX
            derivLeaddbsDir = fileparts(options.subj.subjDir); % .../derivatives/leaddbs
            derivDir = fileparts(derivLeaddbsDir); % .../derivatives
            datasetRoot = fileparts(derivDir); % dataset root
            rawDataDir = fullfile(datasetRoot, 'rawdata', ['sub-', options.subj.subjId]);
            fprintf('Looking in: %s\n', rawDataDir);

            rawAnatFiles = dir(fullfile(rawDataDir, '**', '*_T1w.nii.gz'));
            if isempty(rawAnatFiles)
                rawAnatFiles = dir(fullfile(rawDataDir, '**', '*_T1w.nii'));
            end
            if isempty(rawAnatFiles)
                rawAnatFiles = dir(fullfile(rawDataDir, '**', '*_T2w.nii.gz'));
            end
            if isempty(rawAnatFiles)
                rawAnatFiles = dir(fullfile(rawDataDir, '**', '*_T2w.nii'));
            end

            if ~isempty(rawAnatFiles)
                srcAnat = fullfile(rawAnatFiles(1).folder, rawAnatFiles(1).name);
                [~, anatBaseName, anatExt] = fileparts(rawAnatFiles(1).name);
                if strcmp(anatExt, '.gz')
                    [~, anatBaseName] = fileparts(anatBaseName);
                end
                targetAnat = fullfile(preprocAnatDir, [anatBaseName, '.nii']);

                fprintf('Copying anatomical data from rawdata to preprocessing/anat...\n');
                if endsWith(srcAnat, '.gz')
                    copyfile(srcAnat, [targetAnat, '.gz']);
                    gunzip([targetAnat, '.gz']);
                    delete([targetAnat, '.gz']);
                else
                    copyfile(srcAnat, targetAnat);
                end
                fprintf('Anatomical data copied successfully.\n');

                anatFiles = dir(fullfile(preprocAnatDir, '*_T1w.nii'));
                if isempty(anatFiles)
                    anatFiles = dir(fullfile(preprocAnatDir, '*_T2w.nii'));
                end
            end
        end

        % Filter out SPM tissue class files (c1, c2, c3, etc.) and preprocessed files (r*, sr*, mean*)
        if ~isempty(anatFiles)
            validFiles = {};
            for i = 1:length(anatFiles)
                fname = anatFiles(i).name;
                % Keep only if it doesn't start with c, r, sr, mean
                if ~startsWith(fname, {'c', 'r', 'sr', 'mean'})
                    validFiles{end+1} = fullfile(anatFiles(i).folder, anatFiles(i).name);
                end
            end
            if ~isempty(validFiles)
                options.prefs.prenii_unnormalized = strrep(validFiles{1}, [options.subj.subjDir, filesep], '');
            end
        end

        if ~isfield(options.prefs, 'prenii_unnormalized') || isempty(options.prefs.prenii_unnormalized)
            % No preprocessed anat found - try coregistration/anat directory
            coregAnatDir = fullfile(options.subj.subjDir, 'coregistration', 'anat');
            anatFiles = dir(fullfile(coregAnatDir, '*space-anchorNative*T1w.nii'));
            if isempty(anatFiles)
                anatFiles = dir(fullfile(coregAnatDir, '*space-anchorNative*T2w.nii'));
            end
            if ~isempty(anatFiles)
                anatPath = fullfile(anatFiles(1).folder, anatFiles(1).name);
                options.prefs.prenii_unnormalized = strrep(anatPath, [options.subj.subjDir, filesep], '');
            else
                warning('No anatomical files found in %s or %s - using default', preprocAnatDir, coregAnatDir);
            end
        end
    catch ME
        warning('Failed to set BIDS anatomical path: %s', ME.message);
    end

    % Final check: Ensure b0 path is set correctly (fallback if above logic failed)
    if ~isfield(options.prefs, 'b0') || strcmp(options.prefs.b0, 'b0.nii')
        % BIDS path not set, try to find b0 file
        b0Files = dir(fullfile(options.subj.subjDir, 'preprocessing', 'dwi', '*_b0.nii'));
        if ~isempty(b0Files)
            options.prefs.b0 = fullfile('preprocessing', 'dwi', b0Files(1).name);
        end
    end
else
    error('Not a BIDS dataset!')
end
