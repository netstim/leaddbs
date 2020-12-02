function ea_menu_addsubmit(handles)

s = uimenu('Label','Submit');
cscripts=dir([ea_getearoot,'cluster',filesep,'ea_run_*.m']);
for cscript=1:length(cscripts)
    clustername=ea_underscore2space(cscripts(cscript).name(8:end-2));
    [~,functionname]=fileparts(cscripts(cscript).name);
    c=uimenu(s,'Label',clustername);
    uimenu(c,'Label','Submit Jobs','Callback',{@ea_run_cluster,functionname,handles,'run'});
    uimenu(c,'Label','Export code that submits jobs','Callback',{@ea_run_cluster,functionname,handles,'export'});
end
