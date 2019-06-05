function ea_csremovedrawings(handles)



lineobj=getappdata(handles.checkstructures,'lineobj');
delete(lineobj);
arrhandles=getappdata(handles.checkstructures,'arrhandles');
for arr=1:length(arrhandles)
delete(arrhandles{arr});
end


tp_plots=getappdata(handles.checkstructures,'tp_plots');
for plt=1:length(tp_plots)
delete(tp_plots{plt});
end

