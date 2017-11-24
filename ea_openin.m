function ea_openin(hobj,~,appname,handles)

callingfig=hobj.Parent.Parent.Name;
if strcmp(callingfig,'Lead Group Analysis')
   M=getappdata(handles.leadfigure,'M'); 
    uipatdir=M.patient.list;
    if ~iscell(uipatdir)
        uipatdir={uipatdir};
    end
    if ~isempty(uipatdir)
        uipatdir=uipatdir{get(handles.patientlist,'Value')};
    else
        uipatdir=[];
    end
else
    uipatdir=getappdata(handles.leadfigure,'uipatdir');
end

if ~isempty(uipatdir)
    feval(appname,'loadsubs',uipatdir);
else
    feval(appname);
end

