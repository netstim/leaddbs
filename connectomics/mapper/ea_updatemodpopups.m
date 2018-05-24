function ea_updatemodpopups(mdl,sf,handles)

set(handles.fiberspopup,'String',mdl(sf==1));
set(handles.fmripopup,'String',mdl(sf==2));
if isempty(get(handles.fiberspopup,'String'))
    set(handles.fiberspopup,'String','No structural connectome found.');
end
if isempty(get(handles.fmripopup,'String'))
    set(handles.fmripopup,'String','No functional connectome found.');
end

ea_resetpopup(handles.fmripopup);
