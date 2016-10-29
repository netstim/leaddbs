function ea_resetpopup(pophandle)

if ~iscell(get(pophandle,'String'))
    set(pophandle,'Value',1);
    return
end

if get(pophandle,'Value')>length(get(pophandle,'String'))
    set(pophandle,'Value',length(get(pophandle,'String')));
end