function ea_storeui(handles)

try
    chooseboxname=get(handles.patdir_choosebox,'String');
catch
    return
end

% determine if patientfolder is set
switch chooseboxname
    case 'Choose Patient Directory'
        outdir=ea_getearoot;
    otherwise
        if length(chooseboxname)>=8 && strcmp(chooseboxname(1:8),'Multiple')
        	outdir=ea_getearoot;
        else
            outdir=[get(handles.patdir_choosebox,'String'),filesep];
        end
end

try % only works when calling from core lead (not lead_connectome)
    ea_updatestatus(handles);
end

options=ea_handles2options(handles);
try
    save([outdir,'ea_ui'],'-struct','options');
end

