function ea_openin(hobj,~,appname,handles)


callingfig=hobj.Parent.Parent.Name;
if strcmp(callingfig,'Lead-Group Analysis');
   M=getappdata(handles.leadfigure,'M'); 
    uipatdir=M.patient.list;
    if ~iscell(uipatdir)
        uipatdirc={uipatdir};
        uipatdir=uipatdirc;
    end
    uipatdir=uipatdir{get(handles.patientlist,'Value')};
else
    
    uipatdir=getappdata(handles.leadfigure,'uipatdir');
end

feval(appname,'loadsubs',uipatdir);

