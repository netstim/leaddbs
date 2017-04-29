function  setRelPosition( hObject,units,hObject_pos )
%this function invokes set(hObject,'Position',hObject_pos) while with 'units' provided

old_units=get(hObject,'units');
set(hObject,'units',units);
set(hObject,'Position',hObject_pos);
set(hObject,'units',old_units);

end

