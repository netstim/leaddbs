function ea_load_pts(handles,uipatdir,patsub)

if ~exist('patsub','var')
    patsub='patients';
end

if length(uipatdir)>1
    set(handles.patdir_choosebox,'String',['Multiple (',num2str(length(uipatdir)),')']);
    set(handles.patdir_choosebox,'TooltipString',ea_strjoin(uipatdir,', '));
else
    set(handles.patdir_choosebox,'String',uipatdir{1});
    set(handles.patdir_choosebox,'TooltipString',uipatdir{1});
end

% store patient directories in figure


setappdata(handles.leadfigure,'uipatdir',uipatdir);
try
ea_switchctmr(handles);
end

ea_getui(handles); % update ui from patient
ea_storeui(handles); % save in pt folder
ea_addrecentpatient(handles,uipatdir,['Recent ',patsub,':'],patsub);