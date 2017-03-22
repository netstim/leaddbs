function ea_switchctmr(varargin)
% varargin: handles, switchto (1=MR, 2=CT; optional).
% being called when new patient is loaded.
handles=varargin{1};
switchto=0;

if nargin==1 % autodetect
    % check if MRCT popup is set correctly
    modality=ea_checkctmrpresent(handles);
    switchto=find(modality);
    if any(modality)
        try
        if ~modality(get(handles.MRCT,'Value'))
            set(handles.MRCT,'Value',switchto); 
        end
        end
    end
else
    % switch MRCT popup to specified handle.
    switchto=varargin{2};
end

if ~(sum(switchto>0)>1) && ~isempty(switchto) % e.g. MR and CT present
    switch switchto
        case 1 % MR
                set(handles.coregct_checkbox,'Enable','off');
                set(handles.coregct_checkbox,'Value',0);
                set(handles.coregctmethod,'Enable','off');
                set(handles.coregctcheck,'Enable','off');
                set(handles.coregctcheck,'Value',0);
                set(handles.coregthreshs,'Enable','off');
        case 2 % CT
                set(handles.coregct_checkbox,'Enable','on');
                set(handles.coregctmethod,'Enable','on');
                set(handles.coregctcheck,'Enable','on');
                set(handles.coregthreshs,'Enable','on');
    end
end
ea_updatestatus(handles);