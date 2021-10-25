function ea_load_pts(handles,uipatdir)

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
        if isfile(fullfile(uipatdir{1}, 'ea_ui.mat')) % Old dataset detected
            msg = {'{\bfOld dataset detected, would you like to migrate it to BIDS?}';
                   ['Thank you for your interest in Lead-DBS! Since version 2.6, we have re-organized the way Lead-DBS acesses and stores data.' ,...
                   'This implies changes to the organization of your input and output data. The main objective to set standards for data organization was to promote' ,...
                   'data sharing and open science initiatives. For more information and details on specific changes, please refer to our manual [insert url]. ' ,...
                   'lead-import is a tool developed to automatically assist you in moving your dataset from the classic lead-dbs to the bidsified version. ',...
                   '{\bfIf you wish to run BIDS import tool, please click on ''Yes''. Otherwise, you will not be able to use Lead-DBS.}']};
            opts.Default = 'Cancel';
            opts.Interpreter = 'tex';
            choice = questdlg(msg, '', 'Yes', 'Cancel', opts);
            if strcmp(choice, 'Yes')
                options.prefs = ea_prefs;
                waitfor(lead_import(uipatdir, options, handles));
                BIDSRoot = getappdata(handles.leadfigure,'BIDSRoot');
                subjId = getappdata(handles.leadfigure,'subjID');
                uipatdir = {fullfile(BIDSRoot, 'derivatives', 'leaddbs', ['sub-', subjId{1}])};
            else
                return;
            end
        elseif isfolder(fullfile(uipatdir{1}, 'dicom')) ... % DICOM folder detected
                || isfolder(fullfile(uipatdir{1}, 'DICOM')) ...
                || endsWith(uipatdir{1}, 'dicom', 'IgnoreCase',true)
            msg = {'{\bfDICOM folder found, will run DICOM to NIfTI conversion!}'};
            opts.Interpreter = 'tex';
            opts.WindowStyle = 'modal';
            waitfor(msgbox(msg, '', 'help', opts));
            options.prefs = ea_prefs;
            waitfor(lead_import(uipatdir, options, handles));
            BIDSRoot = getappdata(handles.leadfigure,'BIDSRoot');
            subjId = getappdata(handles.leadfigure,'subjID');
            uipatdir = {fullfile(BIDSRoot, 'derivatives', 'leaddbs', ['sub-', subjId{1}])};
        else
            error('Neither BIDS dataset nor patient folder found!');
        end
    elseif isBIDSRoot % Is BIDS root folder
        BIDSRoot = uipatdir{1};
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
            
            if ~iscell(sourceData)
                sourceData = {sourceData};
            end

            % BIDSRoot is the selected folder
            BIDSRoot = uipatdir{1};
            setappdata(handles.leadfigure, 'BIDSRoot', BIDSRoot);
                   
            subjId = regexp(sourceData, ['(?<=sourcedata\', filesep, 'sub-).*'], 'match', 'once');
            % call lead_migrate
            msg = {'{\bfBIDS dataset with sourcedata found, will run DICOM to NIfTI conversion!}'};
            opts.Interpreter = 'tex';
            opts.WindowStyle = 'modal';
            waitfor(msgbox(msg, '', 'help', opts));
            options.prefs = ea_prefs;

            waitfor(lead_import(sourceData, options, handles));
            uipatdir = strrep(sourceData, 'sourcedata', ['derivatives', filesep, 'leaddbs']);
            
        else
            error('BIDS dataset detected but both sourcedata and rawdata folders are empty!');
        end
    end
else % Multiple patient folders, suppose dataset has already been migrated to BIDS
    BIDSRoot = regexp(uipatdir{1}, ['^.*(?=\', filesep, 'derivatives\', filesep, 'leaddbs)'], 'match', 'once');
    if isempty(BIDSRoot)
        error('Please select patient folders under DATASET/derivatives/leaddbs or migrate the dataset to BIDS format first!');
    else
        subjId = regexp(uipatdir, ['(?<=leaddbs\', filesep, 'sub-).*'], 'match', 'once');
    end
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
    set(handles.patdir_choosebox,'String',uipatdir{1});
    set(handles.patdir_choosebox,'TooltipString',uipatdir{1});
else % Multiple patients mode
    set(handles.patdir_choosebox,'String',['Multiple (',num2str(length(uipatdir)),')']);
    set(handles.patdir_choosebox,'TooltipString',ea_strjoin(uipatdir,'\n'));
end

% Initialize BIDS class
bids = BIDSFetcher(BIDSRoot);

% store patient directories in figure
setappdata(handles.leadfigure, 'uipatdir', uipatdir);
setappdata(handles.leadfigure, 'bids', bids);
setappdata(handles.leadfigure, 'subjId', subjId);

% Set up MR/CT popupmenu and status text
if isfield(handles, 'MRCT')
    ea_switchctmr(handles);
end

if ~ismember(handles.prod, {'mapper'})
    ea_getui(handles); % update ui from patient
end

ea_storeui(handles); % save in pt folder

ea_addrecent(handles, uipatdir, 'patients');

% Return when BIDS dataset is not yet ready
if ~getappdata(handles.leadfigure, 'rawImageJSONExist')
    return;
end

% check if reconstruction is present and assign side-toggles accordingly:
if length(uipatdir) == 1 && isfield(handles, 'side1')
    
    recon = bids.getRecon(subjId{1});
    if isfile(recon.recon)
        load(recon.recon);
        elnum = sum(cellfun(@(f) ~isempty(f), regexp(fieldnames(handles),'^side\d+$','match')));

        % Reset electrode button status
        for el=1:elnum
            set(handles.(['side',num2str(el)]), 'Value', 0);
        end

        % Set electrode button status
        for el=1:length(reco.native.coords_mm)
            if ~isempty(reco.native.markers(el).head)
                set(handles.(['side',num2str(el)]), 'Value', 1);
            end
        end

        try
            elmodel=ea_get_first_notempty_elmodel(reco.props);
            [~,locb] = ismember({elmodel},handles.electrode_model_popup.String);
            set(handles.electrode_model_popup,'Value',locb);
            clear locb
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
        set(handles.seeddefpopup, 'String', [{'Manually choose seeds','Manually choose parcellation'}, commonStims]);
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
