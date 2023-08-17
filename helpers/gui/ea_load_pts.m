function ea_load_pts(handles, uipatdir)

if ~iscell(uipatdir)
    uipatdir = {uipatdir};
end

uipatdir = GetFullPath(uipatdir);

if isLegacyFolder(uipatdir{1})
    % Input folder is legacy patient folder
    if  strcmp(handles.datasetselect.String, 'Choose Dataset Directory')
        options.prefs = ea_prefs;
        BIDSRoot = ea_getdataset(options, handles);
        if ~BIDSRoot
            return;
        end
    end
    BIDSRoot = handles.datasetselect.String;
    subjId = ea_legacy2bids(uipatdir, BIDSRoot, 0);
    if ~iscell(subjId)
        subjId = {subjId};
    end
elseif isBIDSSubjFolder(uipatdir{1})
    % Input folder is BIDS subj folder
    [~, BIDSRoot, subjId] = isBIDSSubjFolder(uipatdir);

    % Patient from another dataset loaded
    if strcmp(handles.prod, 'dbs') && ~ismember(handles.datasetselect.String, {BIDSRoot, 'Choose Dataset Directory'})
        stack = dbstack;
        if any(contains({stack.name}, 'AddPatientButtonPushed'))
            figure(handles.leadfigure);
            answer = uiconfirm(handles.leadfigure,...
                'Selected patients are from another dataset. What would you like to do?',...
                'Add Patient',...
                'Options', {'Copy the selected patients to the current dataset', 'Switch to the selected dataset', 'Cancel'},...
                'DefaultOption', 2, 'CancelOption', 3);
            switch answer
                case 'Copy the selected patients to the current dataset'
                    existedSubjId = subjId(ismember(subjId, handles.patientlist.Data.subjId));
                    if ~isempty(existedSubjId)
                        answer = uiconfirm(handles.leadfigure,...
                            {'These patient IDs already exist in the current dataset:', strjoin(existedSubjId, ', ')}, 'Copying Patients Data...', ...
                             'Options', {'Skip', 'Overwrite', 'Cancel'},...
                            'DefaultOption', 1, 'CancelOption', 3);
                        switch answer
                            case 'Skip'
                                subjId = subjId(~ismember(subjId, handles.patientlist.Data.subjId));
                            case 'Overwrite'
                                ea_cprintf('CmdWinWarnings', 'Deleting existing patient folders...\n');
                                ea_delete([fullfile(handles.datasetselect.String, 'derivatives', 'leaddbs', strcat('sub-', existedSubjId))
                                    fullfile(handles.datasetselect.String, 'rawdata', strcat('sub-', existedSubjId))
                                    fullfile(handles.datasetselect.String, 'sourcedata', strcat('sub-', existedSubjId))]);
                            case 'Cancel'
                                return;
                        end
                    end

                    if ~isempty(subjId)
                        ea_cprintf('CmdWinWarnings', 'Copying patient folders...\n');
                        for i=1:length(subjId)
                            src = {fullfile(BIDSRoot, 'derivatives', 'leaddbs', ['sub-', subjId{i}])
                                fullfile(BIDSRoot, 'rawdata', ['sub-', subjId{i}])
                                fullfile(BIDSRoot, 'sourcedata', ['sub-', subjId{i}])};
                            dst = {fullfile(handles.datasetselect.String, 'derivatives', 'leaddbs', ['sub-', subjId{i}])
                                fullfile(handles.datasetselect.String, 'rawdata', ['sub-', subjId{i}])
                                fullfile(handles.datasetselect.String, 'sourcedata', ['sub-', subjId{i}])};
                            cellfun(@(x, y) isfolder(x) && copyfile(x,y), src, dst, 'Uni', 0);
                        end
                    else
                        if ~isempty(handles.patientlist.Selection)
                            subjId = handles.patientlist.Data.subjId(handles.patientlist.Selection);
                        else
                            subjId = handles.patientlist.Data.subjId(1);
                        end
                    end

                    BIDSRoot = handles.datasetselect.String;
                case 'Cancel'
                    return;
            end
        end
    end
elseif isBIDSFolder(uipatdir{1})
    % Input folder is BIDS dataset root folder, derivatives folder, rawdata folder or sourcedata folder 
    [~, BIDSRoot, subjId] = isBIDSFolder(uipatdir{1});
else % NIfTI or DICOM folder
    % Check if dataset has already been loaded
    if strcmp(handles.prod, 'dbs')
        if  strcmp(handles.datasetselect.String, 'Choose Dataset Directory')
            options.prefs = ea_prefs;
            BIDSRoot = ea_getdataset(options, handles);
            if ~BIDSRoot
                return;
            end
        end

        BIDSRoot = handles.datasetselect.String;
        subjId = cell(length(uipatdir) ,1);
        
        for i = 1:length(uipatdir)
            if isNIfTIFolder(uipatdir{i})
                [~, subjId{i}] = isNIfTIFolder(uipatdir{i});
                subjId{i} = validateSubjId(subjId{i});
                derivativesFolder = fullfile(BIDSRoot, 'derivatives', 'leaddbs', ['sub-', subjId{i}]);
                rawFolder = fullfile(BIDSRoot, 'rawdata', ['sub-', subjId{i}]);
                ea_mkdir({derivativesFolder; rawFolder});
                copyfile(uipatdir{i}, fullfile(rawFolder, 'unsorted'));
            else
                [checkFlag, subjId{i}] = isDICOMFolder(uipatdir{i});
                subjId{i} = validateSubjId(subjId{i});
                if checkFlag
                    derivativesFolder = fullfile(BIDSRoot, 'derivatives', 'leaddbs', ['sub-', subjId{i}]);
                    rawFolder = fullfile(BIDSRoot, 'rawdata', ['sub-', subjId{i}]);
                    sourceFolder = fullfile(BIDSRoot, 'sourcedata', ['sub-', subjId{i}]);
                    ea_mkdir({derivativesFolder; rawFolder; sourceFolder});

                    if endsWith(uipatdir{i}, 'dicom', 'IgnoreCase', true)
                        ea_mkdir(fullfile(sourceFolder, 'DICOM'));
                        ea_cprintf('CmdWinWarnings','Copying DICOM files...\n');
                        copyfile(fullfile(uipatdir{i}, '*'), fullfile(sourceFolder, 'DICOM'));
                    else
                        ea_cprintf('CmdWinWarnings','Copying DICOM files...\n');
                        copyfile(uipatdir{i}, fullfile(sourceFolder, 'DICOM'));
                    end
                    ea_cprintf('CmdWinWarnings','Done.\n');
                else
                    ea_cprintf('CmdWinWarnings', 'Incompatible folder: %s\n', uipatdir{i});
                end
            end
        end
        
        subjId = subjId(~cellfun(@isempty, subjId));
    end
end

% Initialize BIDS class
bids = BIDSFetcher(BIDSRoot);
setappdata(handles.leadfigure, 'bids', bids);
setappdata(handles.leadfigure, 'subjId', subjId);

if ~isempty(bids.subjId)
    uipatdir = fullfile(BIDSRoot, 'derivatives', 'leaddbs', strcat('sub-', subjId));
else
    uipatdir = {'No Patient Selected'};
end

setappdata(handles.leadfigure, 'uipatdir', uipatdir);

if ~isempty(bids.subjId)
    backgroundColor = repmat([1,1,1], length(bids.subjId), 1);
    handles.patientlist.BackgroundColor = backgroundColor;
end

if strcmp(handles.prod, 'dbs')
    handles.datasetselect.String = BIDSRoot;
    handles.datasetselect.TooltipString = BIDSRoot;
    ea_addrecent(handles, {BIDSRoot}, 'datasets');

    if isempty(bids.subjId)
        handles.patientlist.Data = [];
        ea_cprintf('CmdWinWarnings', 'Empty BIDS dataset found!\n');
        return;
    else
        % Set patient listbox
        handles.patientlist.Data = cell2table(bids.subjId, 'VariableNames', {'subjId'});
        handles.patientlist.Selection = find(ismember(bids.subjId', subjId));

        % Check if there are patients not imported yet
        subjDataOverview = bids.subjDataOverview;
        subjNotImported = subjDataOverview.Row(~subjDataOverview.hasRawimagesJson & (subjDataOverview.hasUnsortedRawdata | subjDataOverview.hasSourcedata));
        if ~isempty(subjNotImported)
            highlightInd = find(ismember(bids.subjId', subjNotImported));
            % Highlight only the subjects not imported
            backgroundColor(highlightInd,:) = repmat([249/255,174/255,174/255], length(highlightInd), 1);
            handles.patientlist.BackgroundColor = backgroundColor;

            arrayfun(@(x) set(x, 'Value', 0) , findobj(handles.registrationtab, 'Type', 'uicheckbox'));
            arrayfun(@(x) set(x, 'Value', 0) , findobj(handles.localizationtab, 'Type', 'uicheckbox'));
            arrayfun(@(x) set(x, 'Value', 0) , findobj(handles.optionaltab, 'Type', 'uicheckbox'))

            handles.processtabgroup.SelectedTab = handles.importtab;
            if any(subjDataOverview.hasSourcedata(subjNotImported))
                handles.dicom2bidscheckbox.Value = 1;
            else
                handles.dicom2bidscheckbox.Value = 0;
            end
            if any(subjDataOverview.hasUnsortedRawdata(subjNotImported))
                handles.nifti2bidscheckbox.Value = 1;
            else
                handles.nifti2bidscheckbox.Value = 0;
            end

            handles.statusone.String = 'Unsorted NIfTI/DICOM folder found (highlighted), please run Import first.';
            handles.statustwo.String = '';
        else
            handles.dicom2bidscheckbox.Value = 0;
            handles.nifti2bidscheckbox.Value = 0;
            stack = dbstack;
            if ~ismember('ea_run', {stack.name})
                handles.processtabgroup.SelectedTab = handles.registrationtab;
            end
            handles.statusone.String = '';
            handles.statustwo.String = '';
        end

        % Set Add DICOMs/NIfTIs Menu/Button Visibility
        if length(handles.patientlist.Selection) == 1
            handles.AddDICOMsButton.Visible = 'on';
            handles.AddNIfTIsButton.Visible = 'on';
            handles.UpdateNIfTIsButton.Visible = 'on';
            handles.AddDICOMsMenu.Visible = 'on';
            handles.AddNIfTIsMenu.Visible = 'on';
            handles.UpdateNIfTIsMenu.Visible = 'on';
        else
            handles.AddDICOMsButton.Visible = 'off';
            handles.AddNIfTIsButton.Visible = 'off';
            handles.UpdateNIfTIsButton.Visible = 'off';
            handles.AddDICOMsMenu.Visible = 'off';
            handles.AddNIfTIsMenu.Visible = 'off';
            handles.UpdateNIfTIsMenu.Visible = 'off';
        end

        % Return when there are patients not imported yet
        if ~isempty(subjNotImported)
            return;
        end
    end
end

% Update ui from patient
if ~ismember(handles.prod, {'mapper'})
    ea_getui(handles);
else
    if length(uipatdir) > 1
        handles.patdir_choosebox.String = ['Multiple (', num2str(length(uipatdir)), ')'];
        handles.patdir_choosebox.TooltipString = strjoin(uipatdir, '\n');
    else
        handles.patdir_choosebox.String = uipatdir{1};
        handles.patdir_choosebox.TooltipString = uipatdir{1};
    end
end

% Set up MR/CT popupmenu and status text
if isfield(handles, 'MRCT')
    ea_switchctmr(handles);
end

ea_storeui(handles); % save in pt folder

ea_addrecent(handles, {BIDSRoot}, 'datasets');
ea_addrecent(handles, uipatdir, 'patients');

% check if reconstruction is present and assign side-toggles accordingly:
if length(uipatdir) == 1 && isfield(handles, 'side1')
    if bids.checkModality(subjId{1}, bids.settings.preferMRCT) ~= 3
        recon = bids.getRecon(subjId{1});
        if isfile(recon.recon)
            load(recon.recon);
            elnum = sum(cellfun(@(f) ~isempty(f), regexp(fieldnames(handles),'^side\d+$','match')));
    
            % Reset electrode button status
            for el=1:elnum
                set(handles.(['side',num2str(el)]), 'Value', 0);
            end
    
            % Set electrode button status
            if isfield(reco, 'native')
                recoType = 'native';
            elseif isfield(reco, 'mni')
                recoType = 'mni';
            else
                recoType = '';
            end

            if ~isempty(recoType)
                for el=1:length(reco.(recoType).coords_mm)
                    if ~isempty(reco.(recoType).markers(el).head)
                        set(handles.(['side',num2str(el)]), 'Value', 1);
                    end
                end
            end

            try
                elmodel = ea_get_first_notempty_elmodel(reco.props);
                uiprefsFile = bids.getPrefs(subjId{1}, 'uiprefs', 'mat');
                uiprefs = load(uiprefsFile);
                if ~strcmp(uiprefs.elmodel, elmodel)
                    ea_cprintf('CmdWinWarnings', ...
                        ['Chosen electrode in the GUI (%s) for patient "%s" doesn''t match the stored reconstruction (%s)!\n', ...
                        'Reset to model in the stored recontruction now. Please rerun "Localize DBS electrodes" in case it''s not the correct model.\n'], ...
                        uiprefs.elmodel, subjId{1}, elmodel);

                    try
                        handles.electrode_model_popup.Value = find(ismember(handles.electrode_model_popup.String, elmodel));
                    catch
                         ea_cprintf('CmdWinErrors', ...
                             ['The stored electrode model (%s) for patient "%s" seems not compatible in the current release of LeadDBS.\n', ...
                             'Please double check and choose a correct model.\n'], elmodel, subjId{1});
                    end

                    uiprefs.elmodel = elmodel;
                    save(uiprefsFile, '-struct', 'uiprefs');

                    handles.processtabgroup.SelectedTab = handles.localizationtab;
                    % uialert(handles.leadfigure, ...
                    %     sprintf(['Chosen electrode (%s) for patient "%s" doesn''t match stored reconstruction (%s)! ', ...
                    %     'Please rerun "Localize DBS electrodes".\n'], uiprefs.elmodel, subjId{1}, elmodel), "Warning");
                else
                    [~,locb] = ismember({elmodel},handles.electrode_model_popup.String);
                    set(handles.electrode_model_popup, 'Value', locb);
                    clear locb
                end
            catch
                ea_cprintf('CmdWinWarnings', 'Failed to determine the electrode model for patient "%s"\n', subjId{1});
            end
        end
    end
end

% add VATs to seeds for connectome mapper or predict case
if isfield(handles,'seeddefpopup')
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
    if strcmp(handles.prod, 'mapper')
        commonStims =strcat('Use VAT:', {' '}, commonStims);
        set(handles.seeddefpopup, 'String', [{'Manually choose seeds'; 'Manually choose parcellation'}; commonStims]);
    else
        set(handles.seeddefpopup, 'String', commonStims);
    end
    ea_resetpopup(handles.seeddefpopup);

    % update cons
    if ~strcmp(get(handles.patdir_choosebox,'String'), 'Choose Patient Directory')
        directory = [uipatdir{1}, filesep];
        selectedparc = 'nan';
        options = ea_handles2options(handles);
        options.prefs = bids.settings;
        [mdl,sf] = ea_genmodlist(directory, selectedparc, options);
        ea_updatemodpopups(mdl, sf, handles);
    end
end


function checkFlag = isLegacyFolder(inputFolder)

checkFlag = isfile(fullfile(inputFolder, 'ea_ui.mat')) || isfile(fullfile(inputFolder, 'ea_reconstruction.mat'));


function [checkFlag, BIDSRoot, subjId] = isBIDSFolder(inputFolder)
% Check if input folder is BIDS dataset root folder, derivatives folder, rawdata folder or sourcedata folder 

inputFolder = erase(inputFolder, filesep + textBoundary('end'));

isBIDSSubfolder = endsWith(inputFolder, {'derivatives', ['derivatives', filesep, 'leaddbs'], 'rawdata', 'sourcedata'});
isBIDSRootfolder = ~isempty(ea_regexpdir(inputFolder, '^(derivatives|rawdata|sourcedata)$', 0, 'd'));

checkFlag = isBIDSSubfolder || isBIDSRootfolder;

BIDSRoot = '';
subjId = {''};

if checkFlag
    if nargout > 1
        if isBIDSRootfolder
            BIDSRoot = inputFolder;
        elseif isBIDSSubfolder
            BIDSRoot = regexp(inputFolder, ['.*(?=\',filesep,'(derivatives|derivatives\',filesep,'leaddbs|rawdata|sourcedata))'], 'match', 'once');
        end

        % Select the first subjId by default
        bids = BIDSFetcher(BIDSRoot);
        if ~isempty(bids.subjId)
            subjId = bids.subjId(1);
        end
    end
end


function [checkFlag, BIDSRoot, subjId] = isBIDSSubjFolder(inputFolder)
% Check if input folder is subj folder with in BIDS derivatives folder, rawdata folder or sourcedata folder 

if ischar(inputFolder)
    inputFolder = {inputFolder};
end

inputFolder = erase(inputFolder, filesep + textBoundary('end'));

checkFlag = false;
subjId = {''};

% Extract BIDS root folder. Suppose folders are from the same dataset in case of cell input.
BIDSRoot = regexp(inputFolder{1}, ['.*(?=\', filesep, '(derivatives\', filesep, 'leaddbs|rawdata|sourcedata)\', filesep, 'sub-[^\W_]+$)'], 'match', 'once');

if ~isempty(BIDSRoot)
    checkFlag = true;
    subjId = regexp(inputFolder, '(?<=sub-)[^\W_]+$', 'match', 'once');
end


function [checkFlag, subjId] = isNIfTIFolder(inputFolder)
% Check if input folder is NIfTI folder

inputFolder = erase(inputFolder, filesep + textBoundary('end'));

checkFlag = ~isempty(ea_regexpdir(inputFolder, '.*\.nii(\.gz)?$', 0, 'f'));

[~, subjId] = fileparts(inputFolder);


function [checkFlag, subjId] = isDICOMFolder(inputFolder)
% Check if input folder is NIfTI folder

inputFolder = erase(inputFolder, filesep + textBoundary('end'));

checkFlag = false;
subjId = '';

if ea_dcmquery(inputFolder) > 0
    checkFlag = true;
    if strcmpi(inputFolder, 'dicom')
        % 'DICOM' folder detected, use parent folder name as subjId
        [~, subjId] = fileparts(fileparts(inputFolder));
    else
        % Use folder name as subjId
        [~, subjId] = fileparts(inputFolder);
    end
end

function subjId = validateSubjId(subjId)
if ~isempty(regexp(subjId, '[\W_]', 'once'))
    subjId = regexprep(subjId, '[\W_]', '');
    ea_cprintf('CmdWinWarnings', 'It looks like you have special chars in your subj folder name.\nWe will use a cleaned name ''%s'' for the BIDS dataset. Please check manually.\n', subjId);
end 
