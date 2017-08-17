function ea_run_cluster(~, ~, clusterfunctionname, handles, cmd)

ea_busyaction('on', handles.leadfigure, 'lead');

options = ea_handles2options(handles);
options.macaquemodus = getappdata(handles.leadfigure,'macaquemodus');

% handle 'prefs', 'lc' and 'd2' options
options.prefs = ea_prefs('');
options = ea_amendtoolboxoptions(options);

% mark the lead task as exported job (cluster run is also marked as
% exported job here)
options.exportedJob = 1;

if isfield(handles,'prod')
    if strcmp(handles.prod,'connectome')
        ea_savelcopts(handles);
    end
end

uipatdirs = getappdata(handles.leadfigure,'uipatdir');

switch cmd
    case 'run'
        for pat = 1:length(uipatdirs)   % submit the job one by one
            % set patient specific options
            options.root = [fileparts(uipatdirs{pat}),filesep];
            [~, options.patientname] = fileparts(uipatdirs{pat});
            options.uipatdirs = uipatdirs(pat);
            % run main function
            feval(eval(['@',clusterfunctionname]),options);
        end

    case 'export'
        options.uipatdirs = uipatdirs;
        ea_export(options,clusterfunctionname);
end

ea_busyaction('off', handles.leadfigure, 'lead');
