function ea_run(cmd, options)
% This function is the main execution function of Lead-DBS. It is
% distributed within Lead-DBS toolbox (www.lead-dbs.org)
% __________________________________________________________________________________
% Copyright (C) 2016 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

if strcmp(cmd, 'runcluster')
    load(options);  % 'options' is the job .mat file in the case
    cmd = 'run';    % rewrite 'cmd' to 'run'
    if ~isdeployed
        ea_setpath; % set path (LEAD and SPM) if needed
        addpath(genpath(options.spmdir));
    end
    needtoexit = 1; % on a cluster, need to close Matlab properly after job is done.
else
    needtoexit = 0;
end

% handle 'prefs', 'lc' and 'd2' options if it's not cluster/exported job
if ~isfield(options, 'exportedJob')
    options.prefs = ea_prefs('');
    options = ea_amendtoolboxoptions(options);
end

if options.d3.autoserver && options.d3.write
    choice = questdlg('Are you sure you want to export results to the server automatically?', ...
        'Auto Server-Export', ...
        'Cancel','Yes','Cancel');
    % Handle response
    switch choice
        case 'Cancel'
            return
    end
end

clc
uipatdirs = options.uipatdirs;

if isempty(uipatdirs)
    uipatdirs = {'No Patient Selected'};
end

if ischar(uipatdirs)
    uipatdirs = {uipatdirs};
end

% Check for special characters in the path
ea_checkSpecialChars(uipatdirs);

% do parallel processing if available and set in ea_prefs.
if length(uipatdirs)>1 && ~isempty(which('parpool')) && options.prefs.pp.do && ~strcmp(cmd,'export')
    try
        delete(gcp('nocreate'));
    end

    pp = parpool(options.prefs.pp.profile, options.prefs.pp.csize);

    opts = cell(1, length(uipatdirs));
    for i = 1:length(opts)
        % set patient specific options
        opts{i} = options;
        opts{i}.root = [fileparts(uipatdirs{i}), filesep];
        [~, opts{i}.patientname] = fileparts(uipatdirs{i});
    end

    parfor i = 1:length(uipatdirs)
        % run main function
        try
            ea_autocoord(opts{i});
        catch
            warning([opts{i}.patientname,' failed. Please run this patient again and adjust parameters. Moving on to next patient.']);
        end
    end

    % Refresh GUI after all process is done
    if strcmp(options.leadprod, 'dbs')
        ea_load_pts(getappdata(options.leadfigure, 'handles'), getappdata(options.leadfigure, 'uipatdir'));
    end

    delete(pp);
else
    switch cmd
        case 'export'
            % mark the lead task as exported job
            options.exportedJob = 1;
            ea_export(options);
        case 'run'
            for i = 1:length(uipatdirs)

                % set patient specific options
                options.root = [fileparts(uipatdirs{i}),filesep];
                [~, options.patientname] = fileparts(uipatdirs{i});
                options.subjInd=i;
                % run main function
                if length(uipatdirs) > 1 % multi mode. Dont stop at errors.
                    try
                        ea_autocoord(options);
                    catch
                        warning([options.patientname,' failed. Please run this patient again and adjust parameters. Moving on to next patient.' ]);
                    end
                else
                    ea_autocoord(options);
                end
            end

            % Uncheck CheckBox after successful run
            % handles = getappdata(options.leadfigure, 'handles');
            % checkBoxes = findobj(handles.importtab.Children, 'Type', 'uicheckbox');
            % checkBoxes = [checkBoxes; findobj(handles.registrationtab.Children, 'Type', 'uicheckbox')];
            % checkBoxes = [checkBoxes; findobj(handles.localizationtab.Children, 'Type', 'uicheckbox')];
            % checkBoxes = [checkBoxes; findobj(handles.optionaltab.Children, 'Type', 'uicheckbox')];
            % arrayfun(@(x) set(x, 'Value', 0), checkBoxes);

            % Refresh GUI after all process is done
            if strcmp(options.leadprod, 'dbs') && ~isempty(options.uipatdirs)
                ea_load_pts(getappdata(options.leadfigure, 'handles'), getappdata(options.leadfigure, 'uipatdir'));
            end
    end
end

if needtoexit % close Matlab after job is done.
    exit
end
