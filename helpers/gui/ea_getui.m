function ea_getui(handles)

% determine if patientfolder is set
switch get(handles.patdir_choosebox,'String')
    case {'Choose Patient Directory','Multiple'}
        outdir=[ea_getearoot];
    otherwise
        outdir=get(handles.patdir_choosebox,'String');
end

try
    options=load([outdir,'ea_ui']);
    ea_options2handles(options,handles); % update UI
end
