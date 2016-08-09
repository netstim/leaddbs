function ea_run_cluster(~,~,clusterfunctionname,handles)


leadfig=handles.leadfigure;
ea_busyaction('on',leadfig,'lead');


options=ea_handles2options(handles);
options.macaquemodus=getappdata(handles.leadfigure,'macaquemodus');

options.uipatdirs=getappdata(handles.leadfigure,'uipatdir');

feval(eval(['@',clusterfunctionname]),options);

ea_busyaction('off',leadfig,'lead');
