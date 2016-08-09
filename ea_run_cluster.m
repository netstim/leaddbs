function ea_run_cluster(~,~,clusterfunctionname,handles)

leadfig=handles.leadfigure;
ea_busyaction('on',leadfig,'lead');

options=ea_handles2options(handles);
options.macaquemodus=getappdata(handles.leadfigure,'macaquemodus');

uipatdirs=getappdata(handles.leadfigure,'uipatdir');

for pat=1:length(uipatdirs)
    % set patient specific options
    options.root=[fileparts(uipatdirs{pat}),filesep];
    [~,thispatdir]=fileparts(uipatdirs{pat});
    options.patientname=thispatdir;
    options.uipatdirs=uipatdirs(pat); % only process one patient at a time on a cluster (all is submitted).
    % run main function
    feval(eval(['@',clusterfunctionname]),options);
end

ea_busyaction('off',leadfig,'lead');
