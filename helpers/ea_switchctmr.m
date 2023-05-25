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
    set(handles.MRCT, 'TooltipString', '<html>Post-operative image modality (MR/CT/None) will be automatically detected.<br>In case both MR and CT images are present, CT will be chosen by default.<br>You can change this in your preference file by setting ''prefs.preferMRCT'' (1 for MR and 2 for CT).');

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
        case 'None'
            postopModality = 3;
    end

    % Set status text
    ea_updatestatus(handles, subj);
end

if get(handles.MRCT,'Value') ~= postopModality
	set(handles.MRCT, 'Value', postopModality);
end

if  ~strcmp(handles.prod, 'anatomy')
    arrayfun(@(x) set(x, 'Enable', 'on'), handles.optionaltab.Children);
    set(handles.overwriteapproved, 'Enable', 'on');
    switch postopModality
        case 1 % MR
            arrayfun(@(x) set(x, 'Enable', 'on'), handles.registrationtab.Children);
            set(handles.coregctmethod,'Enable','off');
            set(handles.doreconstruction,'Enable','on');
            set(handles.refinelocalization,'Enable','on');
            set(handles.reconmethod,'String',{'TRAC/CORE (Horn 2015)','Manual', 'Slicer (Manual)'});
            % default TRAC/CORE:
            set(handles.reconmethod,'Enable','on');
            if ismember(ea_getspace,{'Waxholm_Space_Atlas_SD_Rat_Brain','MNI_Macaque'})
                set(handles.reconmethod, 'Value', find(ismember(handles.reconmethod.String, 'Manual')));
            else
                set(handles.reconmethod, 'Value', find(ismember(handles.reconmethod.String, bids.settings.reco.method.MRI))); % set to TRAC/CORE algorithm.
            end
            set(handles.targetpopup,'Enable','on');
            set(handles.maskwindow_txt,'Enable','on');
        case 2 % CT
            arrayfun(@(x) set(x, 'Enable', 'on'), handles.registrationtab.Children);
            set(handles.doreconstruction,'Enable','on');
            set(handles.refinelocalization,'Enable','on');
            set(handles.reconmethod,'String',{'Refined TRAC/CORE','TRAC/CORE (Horn 2015)','PaCER (Husch 2017)','Manual', 'Slicer (Manual)'});
            % default PaCER:
            set(handles.reconmethod,'Enable','on');
            if ismember(ea_getspace,{'Waxholm_Space_Atlas_SD_Rat_Brain','MNI_Macaque'})
                set(handles.reconmethod, 'Value', find(ismember(handles.reconmethod.String, 'Manual')));
            else
                set(handles.reconmethod, 'Value', find(ismember(handles.reconmethod.String, bids.settings.reco.method.CT))); % set to PaCER algorithm.
            end
            set(handles.targetpopup,'Enable','off');
            set(handles.maskwindow_txt,'Enable','off');
        case 3 % None
            arrayfun(@(x) set(x, 'Enable', 'on'), handles.registrationtab.Children);
            set(handles.coregctmethod,'Enable','off');
            set(handles.scrf,'Enable','off');
            set(handles.scrf,'Value',0);
            set(handles.doreconstruction,'Enable','off');
            set(handles.refinelocalization,'Enable','off');
            set(handles.reconmethod,'Enable','off');
            set(handles.targetpopup,'Enable','off');
            set(handles.maskwindow_txt,'Enable','off');
    end
end
