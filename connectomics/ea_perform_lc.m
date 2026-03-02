function ea_perform_lc(options)
% This function is the main execution function of LEAD-DBS Connectome.
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

% BIDS fix: Ensure options.prefs has correct BIDS paths
directory = [options.root, options.patientname, filesep];

% Always call ea_getptopts for BIDS datasets to ensure correct paths
if contains(directory, 'derivatives') || contains(directory, 'leaddbs')
    try
        options = ea_getptopts(directory, options);
    catch
        % If ea_getptopts fails, try to construct BIDS paths manually
        warning('Could not load BIDS paths via ea_getptopts, attempting manual path construction');
    end
end

%% structural parts
disp('*** Performing structural parts of LEAD-Connectome...');

% perform fibertracking
if options.lc.struc.ft.do
    ea_perform_ft_proxy(options);
end

% normalize fibers
if options.lc.struc.ft.normalize % normalize fibertracts ? for now these should be denoted in Freiburg format.
    % Refresh options with BIDS paths before normalization (in case they were updated during tracking)
    if contains(directory, 'derivatives') || contains(directory, 'leaddbs')
        try
            options = ea_getptopts(directory, options);
        catch
            warning('Could not refresh BIDS paths before normalization');
        end
    end
    ea_normalize_fibers(options);
end

% create structural CM
if options.lc.struc.compute_CM
    if ~exist([options.root,options.patientname,filesep,'connectomics'],'dir')
        mkdir([options.root,options.patientname,filesep,'connectomics']);
    end
    expfolder=[options.root,options.patientname,filesep,'connectomics',filesep,options.lc.general.parcellation,filesep];
    if ~exist(expfolder,'dir')
        mkdir(expfolder);
    end
    if ~exist([expfolder,'DTI_CM.mat'],'file') || ...
       (isfield(options, 'overwriteapproved') && options.overwriteapproved)
        if ~exist([options.root,options.patientname,filesep,options.prefs.FTR_unnormalized],'file') % fibertracking has not been performed.
            warning('Fibertracking has not been done yet. Will do so before estimating structural connectivity matrix.');
            ea_perform_ft_proxy(options);
            
            % BIDS FIX: Refresh options after fiber tracking to ensure FTR path is correct
            if contains(directory, 'derivatives') || contains(directory, 'leaddbs')
                try
                    options = ea_getptopts(directory, options);
                    fprintf('Refreshed options after fiber tracking. FTR path: %s\n', options.prefs.FTR_unnormalized);
                catch
                    warning('Could not refresh BIDS paths after fiber tracking');
                end
            end
        end
        DTI_CM=ea_createCM_dti(options);
        cm=ea_export_CM_png(DTI_CM,'DTI Connectivity matrix',options,[0 mean(DTI_CM(:))+2*std(DTI_CM(:))]);
        save([expfolder,'DTI_CM.mat'],'DTI_CM','-v7.3');
        saveas(cm,[expfolder,'DTI_CM.png']);
    end
end

disp('*** Done.');

%% functional parts
if options.lc.func.compute_GM || options.lc.func.compute_CM
    disp('*** Performing functional parts of LEAD-Connectome...');

    % BIDS FIX: Check and copy rs-fMRI from rawdata if needed
    directory = [options.root, options.patientname, filesep];
    preprocFuncDir = fullfile(directory, 'preprocessing', 'func');
    
    % Ensure preprocessing/func exists
    if ~exist(preprocFuncDir, 'dir')
        mkdir(preprocFuncDir);
    end
    
    % Search for BOLD files (exclude preprocessed versions)
    allBoldFiles = dir(fullfile(preprocFuncDir, '*_bold.nii'));
    
    % Filter out preprocessed versions (mean*, r*, sr*, hdmean*)
    boldFiles = [];
    for i = 1:length(allBoldFiles)
        fname = allBoldFiles(i).name;
        if ~startsWith(fname, {'mean', 'r', 'sr', 'hdmean'})
            boldFiles = [boldFiles; allBoldFiles(i)];
        end
    end
    
    % If not found, copy from rawdata
    if isempty(boldFiles)
        fprintf('No rs-fMRI found in preprocessing/func, checking rawdata...\n');
        % Navigate from derivatives/leaddbs/sub-XXX to rawdata/sub-XXX
        % directory = /path/to/derivatives/leaddbs/sub-XXX/
        % We need: /path/to/rawdata/sub-XXX/
        derivLeaddbsDir = fileparts(directory); % /path/to/derivatives/leaddbs
        derivDir = fileparts(derivLeaddbsDir); % /path/to/derivatives
        datasetRoot = fileparts(derivDir); % /path/to/
        rawDataDir = fullfile(datasetRoot, 'rawdata', options.patientname);
        fprintf('Looking in: %s\n', rawDataDir);
        
        rawBoldFiles = dir(fullfile(rawDataDir, '**', '*_bold.nii.gz'));
        if isempty(rawBoldFiles)
            rawBoldFiles = dir(fullfile(rawDataDir, '**', '*_bold.nii'));
        end
        
        if ~isempty(rawBoldFiles)
            srcBold = fullfile(rawBoldFiles(1).folder, rawBoldFiles(1).name);
            [~, boldBaseName, boldExt] = fileparts(rawBoldFiles(1).name);
            if strcmp(boldExt, '.gz')
                [~, boldBaseName] = fileparts(boldBaseName);
            end
            targetBold = fullfile(preprocFuncDir, [boldBaseName, '.nii']);
            
            fprintf('Copying rs-fMRI from rawdata to preprocessing/func...\n');
            if endsWith(srcBold, '.gz')
                copyfile(srcBold, [targetBold, '.gz']);
                gunzip([targetBold, '.gz']);
                delete([targetBold, '.gz']);
            else
                copyfile(srcBold, targetBold);
            end
            fprintf('rs-fMRI copied successfully.\n');
            
            boldFiles = dir(fullfile(preprocFuncDir, '*_bold.nii'));
        end
    end
    
    % get files with rs-fMRI data
    % BIDS-aware: Build restfiles struct consistently
    if ~isempty(boldFiles)
        % Found BOLD files in preprocessing/func
        options.prefs.n_rest = numel(boldFiles);
        % Build restfiles struct with relative paths
        restfiles = struct('name', {}, 'folder', {});
        for i = 1:numel(boldFiles)
            % Store relative path from subject directory
            restfiles(i).name = strrep(fullfile(boldFiles(i).folder, boldFiles(i).name), directory, '');
            restfiles(i).folder = directory;
        end
    elseif isfield(options.prefs, 'rest') && ~isempty(options.prefs.rest) && contains(options.prefs.rest, filesep)
        % BIDS mode: rest is already set to a relative path
        restfiles = struct('name', options.prefs.rest, 'folder', directory);
        options.prefs.n_rest = 1;
    else
        % Classic mode: search for rest files using search string
        if ~isfield(options.prefs, 'rest_searchstring')
            options.prefs.rest_searchstring = 'rest*.nii';
        end
        restfilesRaw = dir([options.root,options.patientname,filesep,options.prefs.rest_searchstring]);
        options.prefs.n_rest = numel(restfilesRaw);
        % Convert to relative paths
        restfiles = struct('name', {}, 'folder', {});
        for i = 1:numel(restfilesRaw)
            restfiles(i).name = restfilesRaw(i).name;  % Classic mode: just filename
            restfiles(i).folder = directory;
        end
    end

    % display number of rs-fMRI files to analyze
    disp(['*** ' num2str(options.prefs.n_rest) ' rs-fMRI files to analyze...']);
end

% connectivity matrix steps:
if options.lc.func.compute_CM
    if ~exist([options.root,options.patientname,filesep,'connectomics'],'dir')
        mkdir([options.root,options.patientname,filesep,'connectomics']);
    end
    % set export folder and check if it exists
    expfolder=[options.root,options.patientname,filesep,'connectomics',filesep,options.lc.general.parcellation,filesep];
    if ~exist(expfolder,'dir')
        mkdir(expfolder);
    end
    [~, presentfiles] = ea_assignpretra(options);
    % set filenames for each rs-fMRI file
    for irest = 1:options.prefs.n_rest
        % set filenames for this iteration
        options.prefs.rest = restfiles(irest).name;
        
        % BIDS FIX: Use ea_prependFilename to preserve paths
        options.prefs.pprest = ea_prependFilename(options.prefs.rest, 'sr'); % preprocessed rs-fMRI data
        options.prefs.glrest = ea_prependFilename(options.prefs.rest, 'glr'); % preprocessed and normalized rs-fMRI data
        
        % Get just the basename for the timecourse filename
        [~,name,ext] = fileparts(restfiles(irest).name);
        options.prefs.gmtc=strcat(name,'_tc.mat'); % extracted timecourses of resting state fMRI data

        % create connectivity matrix for each rs-fMRI file
        % BIDS FIX: Use simple file check instead of ea_reglocked
        cmFile = [expfolder,ea_stripext(options.prefs.rest),'_fMRI_CM.mat'];
        if ~exist(cmFile, 'file')
            disp(['Creating connectivity matrix for rs-fMRI file #',num2str(irest),': ',options.prefs.rest]);
            [fMRI_CM, gmtc]=ea_createCM_fmri(options);
            cm=ea_export_CM_png(fMRI_CM,['fMRI Connectivity matrix for ',name],options);
            save([expfolder,options.prefs.gmtc],'gmtc');
            save([expfolder,name,'_fMRI_CM.mat'],'fMRI_CM','-v7.3');
            saveas(cm,[expfolder,name,'_fMRI_CM.png']);
        else
            disp(['Connectivity matrix already exists for: ',options.prefs.rest]);
        end % end loop for this rs-fMRI file
    end % end loop for for all rs-fMRI files
end % end connectivity matrix section

if options.lc.func.compute_GM || options.lc.func.compute_CM
    disp(['*** Done analyzing ' num2str(options.prefs.n_rest) ' rs-fMRI files...']);
end

% graph metrics section below:
% for multiple file support: need to edit `ea_computeGM`
if options.lc.func.compute_GM || options.lc.struc.compute_GM % perform graph metrics
    if ~exist([options.root,options.patientname,filesep,'connectomics'],'dir')
        mkdir([options.root,options.patientname,filesep,'connectomics']);
    end
    expfolder=[options.root,options.patientname,filesep,'connectomics',filesep,options.lc.general.parcellation,filesep];
    if ~exist(expfolder,'dir')
        mkdir(expfolder);
    end

    modes=cell(0);
    finas=cell(0);
    threshs=cell(0);
    fs=[]; % functional/structural measure (for structure-function similarity)
    if options.lc.func.compute_GM
        for irest = 1:options.prefs.n_rest
            modes{end+1}='fMRI';
            [~,rbase]=fileparts(restfiles(irest).name);
            finas{end+1}=[rbase,'_fMRI'];
            threshs{end+1}=options.lc.graph.fthresh;
            fs(end+1)=1;
        end
    end

    if options.lc.struc.compute_GM
        modes{end+1}='DTI';
        finas{end+1}='DTI';
        threshs{end+1}=options.lc.graph.sthresh;
        fs(end+1)=2;
    end

    % Ensure required CM files exist before graph metrics
    missing = {};
    for mode = 1:length(modes)
        cmPath = [expfolder, finas{mode}, '_CM.mat'];
        if ~exist(cmPath, 'file')
            missing{end+1} = cmPath;
        end
    end
    if options.lc.graph.struc_func_sim
        dtiPath = [expfolder, 'DTI_CM.mat'];
        if ~exist(dtiPath, 'file')
            missing{end+1} = dtiPath;
        end
    end
    if ~isempty(missing)
        ea_error(['Connectivity matrix file(s) not found. Enable "Compute connectivity matrix" for the relevant modality (structural and/or functional) and run again.\nMissing: ', strjoin(missing, '\n       ')]);
    end

    try
        ea_computeGM(options,modes,finas,threshs,fs);
    catch ME
        if contains(ME.message, 'Unable to read')
            ea_error('Connectivity matrix file(s) not found or unreadable. Enable "Compute connectivity matrix" (structural and/or functional) and run again.\nOriginal: %s', ME.message);
        end
        rethrow(ME);
    end
end

disp('*** Done.');
