function hObject_pos = getRelPosition( hObject,units )
%this function returns get(hObject,'Position') while with 'units' provided

old_units=get(hObject,'units');
set(hObject,'units',units);
hObject_pos=get(hObject,'Position');
set(hObject,'units',old_units);

end

