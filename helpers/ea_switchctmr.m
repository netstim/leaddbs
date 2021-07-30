function ea_switchctmr(varargin)
% varargin: handles, switchto (1=MR, 2=CT; optional).
% being called when new patient is loaded.
handles=varargin{1};

if nargin==1 % autodetect
    % check if MRCT popup is set correctly
    modality=ea_checkctmrpresent(handles);
    switchto=find(modality);

    if length(switchto) == 2    % both MR and CT image present
        options.prefs = ea_prefs('');
        switchto = options.prefs.preferMRCT;  % set the modality according to 'prefs.preferMRCT'
    end

    if any(modality)
        try
            if get(handles.MRCT,'Value') ~= switchto
                set(handles.MRCT,'Value',switchto);
            end
        end
    end
else
    % switch MRCT popup to specified handle.
    switchto=varargin{2};
end

if ~isempty(switchto) && ~(length(switchto)==2) && ~strcmp(handles.prod, 'anatomy') % e.g., No MR and CT present, both MR and CT present, called from lead_anatomy
    switch switchto
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

try
    ea_updatestatus(handles);
end
