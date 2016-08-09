function ea_run_cluster(~,~,clusterfunctionname,handles)

leadfig=handles.leadfigure;
ea_busyaction('on',leadfig,'lead');

options=ea_handles2options(handles);
options.macaquemodus=getappdata(handles.leadfigure,'macaquemodus');

options.uipatdirs=getappdata(handles.leadfigure,'uipatdir');

for pat=1:length(options.uipatdirs)
    % set patient specific options
    options.root=[fileparts(options.uipatdirs{pat}),filesep];
    [root,thispatdir]=fileparts(options.uipatdirs{pat});
    options.patientname=thispatdir;
    % run main function
    feval(eval(['@',clusterfunctionname]),options);
end

ea_busyaction('off',leadfig,'lead');
