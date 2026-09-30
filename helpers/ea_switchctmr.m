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
            set(handles.reconmethod, 'Value', 1);
            set(handles.reconmethod,'String',{'TRAC/CORE (Horn 2015)','Manual', 'Slicer (Manual)'});
            % Set recon method
            set(handles.reconmethod,'Enable','on');
            if ismember(ea_getspace,{'Waxholm_Space_Atlas_SD_Rat_Brain','MNI_Macaque'})
                set(handles.reconmethod, 'Value', find(ismember(handles.reconmethod.String, 'Manual')));
            else
                if exist('uiprefs', 'var')
                    % Find the index of the matching method (case-insensitive)
                    idx = find(strcmpi(handles.reconmethod.String, uiprefs.reconmethod), 1);
                    if ~isempty(idx)
                        set(handles.reconmethod, 'Value', idx);
                    else
                        set(handles.reconmethod, 'Value', find(ismember(handles.reconmethod.String, bids.settings.reco.method.MRI))); % set to TRAC/CORE algorithm.
                    end
                else
                    set(handles.reconmethod, 'Value', find(ismember(handles.reconmethod.String, bids.settings.reco.method.MRI))); % set to TRAC/CORE algorithm.
                end
            end
            if contains(handles.reconmethod.String{handles.reconmethod.Value}, 'TRAC')
                set(handles.targetpopup,'Enable','on');
                set(handles.maskwindow_txt,'Enable','on');
            else
                set(handles.targetpopup,'Enable','off');
                set(handles.maskwindow_txt,'Enable','off');
            end
        case 2 % CT
            arrayfun(@(x) set(x, 'Enable', 'on'), handles.registrationtab.Children);
            set(handles.doreconstruction,'Enable','on');
            set(handles.reconmethod, 'Value', 1);
            if ~handles.SEEGCheckBox.Value
                set(handles.reconmethod,'String',{'Refined TRAC/CORE','TRAC/CORE (Horn 2015)','PaCER (Husch 2017)','Manual', 'Slicer (Manual)'});
            else
                set(handles.reconmethod,'String',{'LeGUI (Davis 2021)','Manual','Slicer (Manual)'});
            end
            % Set recon method
            set(handles.reconmethod,'Enable','on');
            if ismember(ea_getspace,{'Waxholm_Space_Atlas_SD_Rat_Brain','MNI_Macaque'})
                set(handles.reconmethod, 'Value', find(ismember(handles.reconmethod.String, 'Manual')));
            elseif handles.SEEGCheckBox.Value
                if exist('uiprefs', 'var')
                    % Find the index of the matching method (case-insensitive)
                    idx = find(strcmpi(handles.reconmethod.String, uiprefs.reconmethod), 1);
                    if ~isempty(idx)
                        set(handles.reconmethod, 'Value', idx);
                    else
                        set(handles.reconmethod, 'Value', find(ismember(handles.reconmethod.String, 'LeGUI (Davis 2021)'))); % set to LeGUI algorithm.
                    end
                else
                    set(handles.reconmethod, 'Value', find(ismember(handles.reconmethod.String, 'LeGUI (Davis 2021)'))); % set to LeGUI algorithm.
                end
            else
                if exist('uiprefs', 'var')
                    % Find the index of the matching method (case-insensitive)
                    idx = find(strcmpi(handles.reconmethod.String, uiprefs.reconmethod), 1);
                
                    if ~isempty(idx)
                        set(handles.reconmethod, 'Value', idx);
                    else
                        set(handles.reconmethod, 'Value', find(ismember(handles.reconmethod.String, bids.settings.reco.method.CT))); % set to PaCER algorithm.
                    end
                else
                    set(handles.reconmethod, 'Value', find(ismember(handles.reconmethod.String, bids.settings.reco.method.CT))); % set to PaCER algorithm.
                end
            end
            if contains(handles.reconmethod.String{handles.reconmethod.Value}, 'TRAC')
                set(handles.targetpopup,'Enable','on');
                set(handles.maskwindow_txt,'Enable','on');
            else
                set(handles.targetpopup,'Enable','off');
                set(handles.maskwindow_txt,'Enable','off');
            end
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

    if handles.SEEGCheckBox.Value && strcmp(handles.reconmethod.String{handles.reconmethod.Value}, 'LeGUI (Davis 2021)') || postopModality == 3
        set(handles.electrode_model_popup, 'Enable', 'off');

        for i = 1:15
            set(handles.(['side', num2str(i)]), 'Enable', 'off');
        end

        set(handles.refinelocalization, 'Value', 0);
        set(handles.refinelocalization, 'Enable', 'off');
    else
        set(handles.electrode_model_popup, 'Enable', 'on');

        for i = 1:15
            set(handles.(['side', num2str(i)]), 'Enable', 'on');
        end

        set(handles.refinelocalization, 'Value', 0);
        set(handles.refinelocalization, 'Enable', 'on');
    end
end
