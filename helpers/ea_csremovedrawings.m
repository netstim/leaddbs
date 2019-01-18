function ea_csremovedrawings(handles)



lineobj=getappdata(handles.checkstructures,'lineobj');
delete(lineobj);
arrhandles=getappdata(handles.checkstructures,'arrhandles');
for arr=1:length(arrhandles)
delete(arrhandles{arr});
end


