function ea_subcorticalrefine_menu(~,~,handles)

uipatdirs=getappdata(handles.leadfigure,'uipatdir');

for pt=1:length(uipatdirs)
    
   ea_subcorticalrefine(uipatdirs{pt},handles); 
end
