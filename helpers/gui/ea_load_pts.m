function ea_load_pts(handles,uipatdir,patsub)

if ~exist('patsub','var')
    patsub = 'patients';
end

if ~iscell(uipatdir)
    uipatdir = {uipatdir};
end

isSubjFolder = 0;
isBIDSRoot = 0;

if length(uipatdir) == 1 % Dragged single folder
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
        end
    end

    if ~isSubjFolder && ~isBIDSRoot
        error('Neither BIDS dataset nor patient folder found!');
    elseif isBIDSRoot % Is BIDS root folder
        rawData = ea_regexpdir([uipatdir{1}, filesep, 'rawdata'], 'sub-', 0);
        rawData = regexprep(rawData, ['\', filesep, '$'], '');
        sourceData = ea_regexpdir([uipatdir{1}, filesep, 'sourcedata'], 'sub-', 0);
        sourceData = regexprep(sourceData, ['\', filesep, '$'], '');

        if ~isempty(rawData) % rawdata folder already exists
            uipatdir = strrep(rawData, 'rawdata', ['derivatives', filesep, 'leaddbs']);
            subjId = regexp(rawData, ['(?<=rawdata\', filesep, 'sub-).*'], 'match');
        elseif ~isempty(sourceData) % sourcedata folder exists
            uipatdir = strrep(rawData, 'sourcedata', ['derivatives', filesep, 'leaddbs']);
            subjId = regexp(sourceData, ['(?<=sourcedata\', filesep, 'sub-).*'], 'match');
        else
            error('Both sourcedata and rawdata folders are empty!');
        end
    end
else % Dragged multiple patient folders
    BIDSRoot = regexp(uipatdir{1}, ['^.*(?=\', filesep, 'derivatives\', filesep, 'leaddbs)'], 'match', 'once');
    if isempty(BIDSRoot)
        error('Please select patient folders under DATASET/derivatives/leaddbs!');
    else
        subjId = regexp(uipatdir, ['(?<=leaddbs\', filesep, 'sub-).*'], 'match', 'once');
    end
end

if isBIDSRoot && length(subjId) > 1 % Multiple patients found
    index = listdlg('PromptString','Select Patient', 'ListString', subjId);
    if isempty(index)
        return;
    else
        uipatdir = uipatdir(index);
        subjId = subjId(index);
    end
end

if length(subjId) == 1 % Only one patient selected
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
ea_switchctmr(handles);

ea_getui(handles); % update ui from patient
ea_storeui(handles); % save in pt folder
ea_addrecentpatient(handles, uipatdir, patsub, patsub);

% check if reconstruction is present and assign side-toggles accordingly:
if length(subjId) == 1 && isfield(handles, 'side1')
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
        stims = ea_dir2cell(dir([uipatdir{pt},filesep,'stimulations',filesep,ea_getspace]));
        if ~exist('remstims', 'var')
            commonStims = stims;
        else
            commonStims = intersect(commonStims, stims);
        end
    end

    % for now only check first subject for pt. specific fibers..
    % find out whether mapper or predict were calling
    if strncmp(handles.leadfigure.Name, 'Lead Connectome Mapper', 22)
        commonStims = ea_prependvat(commonStims);
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
        options.prefs = ea_prefs;
        [mdl,sf] = ea_genmodlist(directory, selectedparc, options);
        ea_updatemodpopups(mdl, sf, handles);
    end
end


function remstims=ea_prependvat(remstims)
for rs=1:length(remstims)
    remstims{rs}=['Use VATs: ',remstims{rs}];
end
