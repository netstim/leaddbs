function ea_menu_initmenu(handles,cmd)

menuprobe=getappdata(handles.leadfigure,'menuprobe');
if isempty(menuprobe)
    % tools menu  & edit prefs present in all apps.
    f = uimenu('Label','Tools');
    uimenu(f,'Label','Edit Preferences file...','Callback',{@ea_editprefs});
    if ismember('acpc',cmd)
        uimenu(f,'Label','Convert ACPC/MNI coordinates','Callback',{@ea_acpcquery,handles.leadfigure});
    end
    if ismember('export',cmd)
        e = uimenu(f,'Label','Export');
        uimenu(e,'Label','Export .PDF files for selected patient(s)','Callback',{@ea_exportpat,'PDF',handles});
        uimenu(e,'Label','Export .STL files for selected patient(s)','Callback',{@ea_exportpat,'STL',handles});
        uimenu(e,'Label','Export .PLY files for selected patient(s)','Callback',{@ea_exportpat,'PLY',handles});
    end
    if ismember('cluster',cmd)
        ea_menu_addsubmit(handles);
        setappdata(handles.leadfigure,'menuprobe',1);
    end
    
    % always add install addons
    g = uimenu('Label','Install');
    [list,commands]=ea_checkinstall('list');
    for l=1:length(list)
       uimenu(g,'Label',list{l},'Callback',{@ea_menuinstall,commands{l},0}); 
    end
    
end