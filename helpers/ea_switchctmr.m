function ea_switchctmr(handles, preferMRCT)
% preferMRCT: 1 = MR, 2 = CT

bids = getappdata(handles.leadfigure,'bids');
subjId = getappdata(handles.leadfigure,'subjId');

if length(subjId) > 1 % Mutiple patient mode
    % Disable MR/CT popupmenu
    set(handles.MRCT,'Enable', 'off');
    set(handles.MRCT, 'TooltipString', '<html>Multiple patients are selected.<br>Enable CT to MRI coregistration setting by default.<br>The actual modality will be automatically detected.');

    % Set status text
    statusone = 'Multiple patients are chosen, CT/MR modality will be automatically detected.';
    set(handles.statusone, 'String', statusone);
    set(handles.statusone, 'TooltipString', statusone);

    % This will enable CT coreg setting in multiple patients mode
    postopModality = 2;
else % Only one patient loaded
    % Make sure MR/CT popupmenu is set correctly
    set(handles.MRCT, 'TooltipString', '<html>Post-operative image modality (MR/CT) will be automatically detected.<br>In case both MR and CT images are present, CT will be chosen by default.<br>You can change this in your preference file by setting ''prefs.preferMRCT'' (1 for MR and 2 for CT).');

    % Check MR/CT preference: first check uiprefs, then LeadDBS settings
    if ~exist('preferMRCT', 'var') || isempty(preferMRCT)
        uiprefsFile = bids.getPrefs(subjId{1}, 'uiprefs', 'mat');
        if isfile(uiprefsFile)
            uiprefs = load(uiprefsFile);
            preferMRCT = uiprefs.modality;
        else
            preferMRCT = bids.settings.preferMRCT;
        end
    end

    % Get subj BIDS struct
    subj = bids.getSubj(subjId{1}, preferMRCT);
    if ~subj.rawImageJSONExist
        setappdata(handles.leadfigure, 'rawImageJSONExist', 0);
        return;
    else
        setappdata(handles.leadfigure, 'rawImageJSONExist', 1);
    end

    % Enable MR/CT popupmenu in case both present
    if subj.bothMRCTPresent
        set(handles.MRCT,'Enable', 'on');
    else
        set(handles.MRCT,'Enable', 'off');
    end

    switch subj.postopModality
        case 'MRI'
            postopModality = 1;
        case 'CT'
            postopModality = 2;
    end

    % Set status text
    ea_updatestatus(handles, subj);
end

if get(handles.MRCT,'Value') ~= postopModality
	set(handles.MRCT, 'Value', postopModality);
end

if  ~strcmp(handles.prod, 'anatomy')
    switch postopModality
        case 1 % MR
            set(handles.coregctmethod,'Enable','off');
            set(handles.reconmethod,'String',{'TRAC/CORE (Horn 2015)','Manual', 'Slicer (Manual)'});
            % default TRAC/CORE:
            set(handles.reconmethod,'enable','on');
            if ismember(ea_getspace,{'Waxholm_Space_Atlas_SD_Rat_Brain','MNI_Macaque'})
                set(handles.reconmethod,'Value',2);
            else
                set(handles.reconmethod,'Value',1); % set to TRAC/CORE algorithm.
            end
            set(handles.targetpopup,'enable','on');
            set(handles.maskwindow_txt,'enable','on');
        case 2 % CT
            set(handles.coregctmethod,'Enable','on');
            set(handles.reconmethod,'String',{'Refined TRAC/CORE','TRAC/CORE (Horn 2015)','PaCER (Husch 2017)','Manual', 'Slicer (Manual)'});
            % default PaCER:
            set(handles.reconmethod,'enable','on');
            if ismember(ea_getspace,{'Waxholm_Space_Atlas_SD_Rat_Brain','MNI_Macaque'})
                set(handles.reconmethod,'Value',4);
            else
                set(handles.reconmethod,'Value',3); % set to PaCER algorithm.
            end
            set(handles.targetpopup,'enable','off');
            set(handles.maskwindow_txt,'enable','off');
    end
end
