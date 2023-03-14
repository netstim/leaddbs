function ea_load_pts_dataset(app,uipatdir)
% function branched off of ea_load_pts (used by GUIDE apps) to be used by
% AppDesigner Apps

if ~iscell(uipatdir)
    uipatdir = {uipatdir};
end

uipatdir = GetFullPath(uipatdir);

isSubjFolder = 0;
isBIDSRoot = 0;

if length(uipatdir) == 1 % Single folder
    if contains(uipatdir{1}, ['derivatives', filesep, 'leaddbs']) % Is patient folder under derivatives
        isSubjFolder = 1;
        BIDSRoot = regexp(uipatdir{1}, ['^.*(?=\', filesep, 'derivatives)'], 'match', 'once');
        subjId = regexp(uipatdir{1}, ['(?<=leaddbs\', filesep, 'sub-).*'], 'match');
    elseif contains(uipatdir{1}, {['rawdata', filesep, 'sub-'], ['sourcedata', filesep, 'sub-']}) % rawdata folder has been selected
        isSubjFolder = 1;
        BIDSRoot = regexp(uipatdir{1}, ['^.*(?=\', filesep, '(rawdata|sourcedata))'], 'match', 'once');
        subjId = regexp(uipatdir{1}, ['(?<=(rawdata|sourcedata)\', filesep, 'sub-).*'], 'match');

        subjDerivativesFolder = fullfile(BIDSRoot, 'derivatives', 'leaddbs', ['sub-', subjId{1}]);
        subjRawdataFolder = fullfile(BIDSRoot, 'rawdata', ['sub-', subjId{1}]);
        subjSourcedataFolder = fullfile(BIDSRoot, 'sourcedata', ['sub-', subjId{1}]);

        % now check if derivatives/sub-xx exists
        if isfile(fullfile(subjDerivativesFolder, 'prefs', ['sub-', subjId{1}, '_desc-rawimages.json']))
            ea_cprintf('CmdWinWarnings', 'rawimages.json detected for "sub-%s":\n', subjId{1});
            uipatdir = {fullfile(BIDSRoot, 'derivatives', 'leaddbs', ['sub-', subjId{1}])};
        elseif ~isempty(ea_regexpdir(subjRawdataFolder, '.*\.nii(\.gz)$', 1, 'f')) % rawimages.json missing, but images exist in rawdata folder
            ea_cprintf('CmdWinWarnings', 'rawdata exists but rawimages.json is not present!\n');
            ea_genrawimagesjson(BIDSRoot, subjId{1});
            uipatdir = {subjDerivativesFolder};
        elseif isfolder(subjSourcedataFolder) && ~isempty(ea_regexpdir(subjSourcedataFolder, '.*', 1, 'f'))
            ea_cprintf('CmdWinErrors', 'rawdata folder is empty! Will try to import data from sourcedata folder...\n');
            options.prefs = ea_prefs;
            waitfor(lead_import(subjSourcedataFolder, options, app, BIDSRoot));
            uipatdir = {subjDerivativesFolder};
        end

    else % Check if it's BIDS root folder
        folders = dir(uipatdir{1});
        folders = {folders.name};
        if ismember('sourcedata', folders) || ismember('rawdata', folders)
            isBIDSRoot = 1;
            BIDSRoot = uipatdir{1};
            if ~ismember(folders,'dataset_description.json')
                disp("could not find dataset description file, generating one now...");
                ea_generate_datasetDescription(uipatdir{1}, 'root_folder');
            end
        end
    end

    if ~isSubjFolder && ~isBIDSRoot
        % try to find out what kind of folder structure was passed
        folder_list = dir_without_dots(uipatdir{1});         % do a listing of the immediate directory first
        dcm_in_subfolder_list = dir(fullfile(uipatdir{1}, '**/*.dcm'));     % find out any .dcm files in any subdir
        raw_nifti_in_subfolder_list = {folder_list.name};
        raw_nifti_filter = {'.nii', '.nii.gz'};
        % Old dataset detected
        if isfile(fullfile(uipatdir{1}, 'ea_ui.mat'))
            folder_type = 'legacy_patient_folder';

        % DICOM folder detected
        elseif endsWith(uipatdir{1}, 'dicom', 'IgnoreCase',true)
            uipatdir = {fileparts(uipatdir{1})};
            folder_type = 'patient_folder_dicom_folder';

        % DICOM folder detected inside the folder
        elseif ~isempty(ea_regexpdir(uipatdir{1}, 'DICOM', 0, 'd'))
            folder_type = 'patient_folder_dicom_folder';

        % if DICOMDIR file inside, assume dicoms are present in one of the folders on this level
        elseif any(strcmp('DICOMDIR', {folder_list.name}))
            folder_list(~[folder_list.isdir]' | startsWith({folder_list.name}', {'.', '..'})) = [];
            from = fullfile(uipatdir{1}, {folder_list.name}');
            to = fullfile(uipatdir{1}, 'DICOM');
            ea_mkdir(to);
            parfor f = 1:length(from)
                % if .dcm file is inside subfolder, just move the folder
                movefile(from{f}, to);
            end
            folder_type = 'patient_folder_dicom_folder';

        % .dcm files are found inside one of the subfolders
        elseif ~isempty(dcm_in_subfolder_list)  %
            from = fullfile({dcm_in_subfolder_list.folder}', {dcm_in_subfolder_list.name}');
            to = strrep({dcm_in_subfolder_list.folder}', uipatdir{1}, fullfile(uipatdir{1}, 'DICOM'));
            ea_mkdir(to);
            parfor f=1:length(dcm_in_subfolder_list)
                movefile(from{f}, to{f});
            end

            folder_type = 'patient_folder_dicom_folder';

        % does not have ea_ui.mat, only has niftis
        elseif ~isempty(raw_nifti_in_subfolder_list) && all(endsWith(raw_nifti_in_subfolder_list, raw_nifti_filter))
            folder_type = 'patient_folder_raw_nifti';

        % Otherwise, look for DICOMs
        else
            warning('off', 'backtrace');
            warning('trying to load as DICOM folder...');
            warning('on', 'backtrace');
            movefile(uipatdir{1}, [uipatdir{1}, '_ORIG']);
            mkdir(fullfile(uipatdir{1}, 'DICOM'));
            movefile([uipatdir{1}, '_ORIG'], fullfile(uipatdir{1}, 'ORIG'));

            if ispc
                dcm2niix = fullfile(ea_getearoot, 'ext_libs', 'dcm2nii', 'dcm2niix.exe');
            else
                dcm2niix = fullfile(ea_getearoot, 'ext_libs', 'dcm2nii', ['dcm2niix.', computer('arch')]);
            end

            cmd = [dcm2niix, ' -r y', ' -o ', ea_path_helper(fullfile(uipatdir{1}, 'DICOM')), ' ', ea_path_helper(fullfile(uipatdir{1}, 'ORIG'))];

            % Search for DICOM files and rename to *.dcm
            if ~ispc
                [~, cmdout] = system(['bash -c "', cmd, '"']);
            else
                [~, cmdout] = system(cmd);
            end

            numDICOMs = regexp(cmdout, '(?<=Converted )\d+(?= DICOMs)', 'match', 'once');
            if strcmp(numDICOMs, '0')
                warning('off', 'backtrace');
                warning('%s DICOMs found!', numDICOMs);
                warning('on', 'backtrace');
                folder_type = '';
            else
                fprintf('%s DICOMs found!\n', numDICOMs);
                ea_delete(fullfile(uipatdir{1}, 'ORIG'));
                folder_type = 'patient_folder_dicom_folder';
            end
        end
       
        
        switch folder_type
            case 'legacy_patient_folder'
                options.prefs = ea_prefs;
                msg = sprintf('Old dataset with legacy files detected,\n would you like to migrate it to BIDS?');
                waitfor(ea_selectdataset(msg,app.leadfigure));
                dest_folder = getappdata(app.leadfigure, 'BIDSRoot');
                if ~isempty(dest_folder)
                    if options.prefs.migrate.interactive
                        waitfor(lead_import(uipatdir, options, app,dest_folder));
                    else
                        ea_lead_import(uipatdir,options,app,dest_folder)
                    end
                    BIDSRoot = getappdata(app.leadfigure,'BIDSRoot');
                    subjId = getappdata(app.leadfigure,'subjID');
                    if ~isempty(BIDSRoot) && ~isempty(subjId)
                        uipatdir = {fullfile(BIDSRoot, 'derivatives', 'leaddbs', ['sub-', subjId{1}])};
                    else
                        return
                    end
                else %user pressed cancel in ea_selectdatasets
                    return;
                end
            case  'patient_folder_dicom_folder'
                options.prefs = ea_prefs;
                msg = sprintf('DICOM folder found,\n should we run DICOM to NIfTI conversion?');
                waitfor(ea_selectdataset(msg,app.leadfigure));
                dest_folder = getappdata(app.leadfigure, 'BIDSRoot');
                if ~isempty(dest_folder)
                    if options.prefs.migrate.interactive
                        waitfor(lead_import(uipatdir, options, app, dest_folder));
                    else
                        ea_lead_import(uipatdir,options,app,dest_folder);
                    end
                    BIDSRoot = getappdata(app.leadfigure,'BIDSRoot');
                    subjId = getappdata(app.leadfigure,'subjID');
                    if ~isempty(BIDSRoot) && ~isempty(subjId)
                        uipatdir = {fullfile(BIDSRoot, 'derivatives', 'leaddbs', ['sub-', subjId{1}])};
                    else
                        return
                    end
                else %user pressed cancel in ea_selectdatasets
                    return
                end
            case 'patient_folder_raw_nifti'
                options.prefs = ea_prefs;
                msg = sprintf('Raw dataset with Nifti files [only] detected,\n would you like to migrate it to BIDS?');
                waitfor(ea_selectdataset(msg,app.leadfigure));
                dest_folder = getappdata(app.leadfigure, 'BIDSRoot');
                if ~isempty(dest_folder)
                    if options.prefs.migrate.interactive
                        waitfor(lead_import(uipatdir, options, app));
                    else
                        dest_folder = ea_uigetdir('*','Select BIDS Dataset folder');
                        if isempty(dest_folder) %user pressed cancel
                            return
                        end
                        ea_lead_import(uipatdir,options,app,dest_folder);
                    end

                    BIDSRoot = getappdata(app.leadfigure,'BIDSRoot');
                    subjId = getappdata(app.leadfigure,'subjID');
                    if ~isempty(BIDSRoot) && ~isempty(subjId)
                        uipatdir = {fullfile(BIDSRoot, 'derivatives', 'leaddbs', ['sub-', subjId{1}])};
                    else
                        return
                    end
                else %user pressed cancel in ea_selectdatasets
                    return;
                end
            otherwise
                error('No compatible files/folders (BIDS dataset/NIfTIs/DICOMs) found!');
        end
    elseif isBIDSRoot % Is BIDS root folder
        BIDSRoot = uipatdir{1};
        ea_checkSpecialChars(BIDSRoot);
        rawData = ea_regexpdir([uipatdir{1}, filesep, 'rawdata'], 'sub-', 0, 'dir');
        rawData = regexprep(rawData, ['\', filesep, '$'], '');
        sourceData = ea_regexpdir([uipatdir{1}, filesep, 'sourcedata'], 'sub-', 0, 'dir');
        sourceData = regexprep(sourceData, ['\', filesep, '$'], '');

        if ~isempty(rawData) % rawdata folder already exists
            uipatdir = strrep(rawData, 'rawdata', ['derivatives', filesep, 'leaddbs']);
            subjId = regexp(rawData, ['(?<=rawdata\', filesep, 'sub-).*'], 'match', 'once');
        elseif ~isempty(sourceData) % sourcedata folder exists
            % in the case of a BIDS dataset root folder as input and sourcedata available for one or more patients
            % trigger DICOM->nii conversion

            % BIDSRoot is the selected folder
            BIDSRoot = uipatdir{1};
            setappdata(app.leadfigure, 'BIDSRoot', BIDSRoot);

            subjId = regexp(sourceData, ['(?<=sourcedata\', filesep, 'sub-).*'], 'match', 'once');
            % call lead_migrate
            msg = {'{\bfBIDS dataset with sourcedata found, will run DICOM to NIfTI conversion!}'};
            opts.Interpreter = 'tex';
            opts.WindowStyle = 'modal';
            waitfor(msgbox(msg, '', 'help', opts));
            options.prefs = ea_prefs;
            dest_folder = BIDSRoot; %source and dest are same
            if options.prefs.migrate.interactive
                waitfor(lead_import(sourceData, options, app));
            else
                ea_lead_import(uipatdir,options,app,dest_folder)
            end
            %
            uipatdir = strrep(sourceData, 'sourcedata', ['derivatives', filesep, 'leaddbs']);

        else
            error('BIDS dataset detected but both sourcedata and rawdata folders are empty!');
        end
    end
else % Multiple patient folders, suppose dataset has already been migrated to BIDS
    BIDSRoot = regexp(uipatdir{1}, ['^.*(?=\', filesep, 'derivatives\', filesep, 'leaddbs)'], 'match', 'once');
    if isempty(BIDSRoot)
        options.prefs = ea_prefs;
        msg = sprintf('Multiple datasets detected,\n would you like to migrate it to BIDS?');
        waitfor(ea_selectdataset(msg,app.leadfigure));
        dest_folder = getappdata(app.leadfigure, 'BIDSRoot');
        if ~isempty(dest_folder)
            if options.prefs.migrate.interactive
                waitfor(lead_import(uipatdir, options, app,dest_folder))
            else
                ea_lead_import(uipatdir,options,app,dest_folder)
            end
            BIDSRoot = getappdata(app.leadfigure,'BIDSRoot');
            subjId = getappdata(app.leadfigure,'subjID');
            derivatives_folder = fullfile(BIDSRoot,'derivatives');
            sourcedata_folder = fullfile(BIDSRoot,'sourcedata');
            rawdata_folder = fullfile(BIDSRoot,'rawdata');
            if isfolder(derivatives_folder)
                uipatdir = fullfile(BIDSRoot,'derivatives','leaddbs',subjId);
            elseif isfolder(sourcedata_folder)
                uipatdir = fullfile(BIDSRoot,'sourcedata',subjId);
            elseif isfolder(rawdata_folder)
                uipatdir = fullfile(BIDSRoot,'rawdata',subjId);

            else
                warning('Something went wrong! Consider migrating the dataset using lead import alone');
                return
            end
            isBIDSRoot = 1;
        else %user pressed cancel in ea_selectdatasets
            return
        end
        %subjId = regexp(uipatdir, ['(?<=leaddbs\', filesep, 'sub-).*'], 'match', 'once');
         %error('Please select patient folders under DATASET/derivatives/leaddbs or migrate the dataset to BIDS format first!');
    else
        subjId = regexp(uipatdir, ['(?<=leaddbs\', filesep, 'sub-).*'], 'match', 'once');
    end

end


if ~iscell(uipatdir)
    uipatdir = {uipatdir};
end

if ~iscell(subjId)
    subjId = {subjId};
end

if isBIDSRoot && length(uipatdir) > 1 % Multiple patients found
    index = listdlg('PromptString','Select Patient', 'ListString', subjId);
    if isempty(index)
        return;
    else
        uipatdir = uipatdir(index);
        subjId = subjId(index);
    end
end

if length(uipatdir) == 1 % Only one patient selected
    set(app.patdir_choosebox,'String',uipatdir{1});
    set(app.patdir_choosebox,'TooltipString',uipatdir{1});
else % Multiple patients mode
    set(app.patdir_choosebox,'String',['Multiple (',num2str(length(uipatdir)),')']);
    set(app.patdir_choosebox,'TooltipString',ea_strjoin(uipatdir,'\n'));
end

% Initialize BIDS class
bids = BIDSFetcher(BIDSRoot);

if ~any(ismember(subjId, bids.subjId))
    % Return when import for all subjs cancelled or failed
    ea_cprintf('CmdWinWarnings', 'Import for sub-%s didn''t go through!\n', subjId{:})
    return;
elseif ~all(ismember(subjId, bids.subjId))
    % Warn if import for some subjs cancelled or failed
    warnSubjId = subjId(~ismember(subjId, bids.subjId));
    ea_cprintf('CmdWinWarnings', 'Import for sub-%s didn''t go through!\n', warnSubjId{:})
end

% store patient directories in figure
setappdata(app.leadfigure, 'uipatdir', uipatdir);
setappdata(app.leadfigure, 'bids', bids);
setappdata(app.leadfigure, 'subjId', subjId);

% Set up MR/CT popupmenu and status text
if isfield(app, 'MRCT')
    ea_switchctmr(app);
end

if ~ismember(app.prod, {'mapper'})
    ea_getui(app); % update ui from patient
end

ea_storeui(app); % save in pt folder

ea_addrecent(app, uipatdir, 'patients');

% check if reconstruction is present and assign side-toggles accordingly:
if length(uipatdir) == 1 && isfield(app, 'side1')
    if bids.checkModality(subjId{1}, bids.settings.preferMRCT) ~= 3
        recon = bids.getRecon(subjId{1});
        if isfile(recon.recon)
            load(recon.recon);
            elnum = sum(cellfun(@(f) ~isempty(f), regexp(fieldnames(app),'^side\d+$','match')));
    
            % Reset electrode button status
            for el=1:elnum
                set(app.(['side',num2str(el)]), 'Value', 0);
            end
    
            % Set electrode button status
            for el=1:length(reco.native.coords_mm)
                if ~isempty(reco.native.markers(el).head)
                    set(app.(['side',num2str(el)]), 'Value', 1);
                end
            end
    
            try
                elmodel=ea_get_first_notempty_elmodel(reco.props);
                [~,locb] = ismember({elmodel},app.electrode_model_popup.String);
                set(app.electrode_model_popup,'Value',locb);
                clear locb
            end
        end
    end
end

% add VATs to seeds for connectome mapper or predict case
if isfield(app,'seeddefpopup')
    for pt=1:length(uipatdir)
        [~, stims] = fileparts(ea_regexpdir(fullfile(uipatdir{pt}, 'stimulations', ea_getspace), '.*', 0, 'dir'));
        if ischar(stims)
            stims = {stims};
        end

        if ~exist('commonStims', 'var')
            commonStims = stims;
        else
            commonStims = intersect(commonStims, stims);
        end
    end

    % for now only check first subject for pt. specific fibers..
    % find out whether mapper or predict were calling
    if strcmp(app.prod, 'mapper')
        commonStims =strcat('Use VAT:', {' '}, commonStims);
        set(app.seeddefpopup, 'String', [{'Manually choose seeds'; 'Manually choose parcellation'}; commonStims]);
    else
        set(app.seeddefpopup, 'String', commonStims);
    end
    ea_resetpopup(app.seeddefpopup);

    % update cons
    if ~strcmp(get(app.patdir_choosebox,'String'), 'Choose Patient Directory')
        directory = [uipatdir{1}, filesep];
        selectedparc = 'nan';
        options = ea_handles2options(app);
        options.prefs = bids.settings;
        [mdl,sf] = ea_genmodlist(directory, selectedparc, options);
        ea_updatemodpopups(mdl, sf, app);
    end
end
