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
if isunix
    special_characters = cellfun(@(x) regexpi(x,'[^a-z0-9\/_-]*','match'), uipatdirs, 'uni', 0);
else
    special_characters = cellfun(@(x) regexpi(x(3:end),'[^a-z0-9\\_-]*','match'), uipatdirs, 'uni', 0);
end
warntxt = '';
for i = find(~cellfun(@isempty, special_characters))'
    warntxt = [warntxt sprintf('The folder: %s\ncontains the unsopported charaters: ''%s''.\n\n', uipatdirs{i}, strjoin(special_characters{i}, ''','''))];
end
if ~isempty(warntxt)
    warndlg([warntxt 'This can make some routines like normalization fail. Please rename the folders and select them again.'], 'Unsupported characters');
end

% do parallel processing if available and set in ea_prefs.
if length(uipatdirs)>1 && ~isempty(which('parpool')) && options.prefs.pp.do && ~strcmp(cmd,'export')
    try
        delete(gcp('nocreate'));
    end

    pp = parpool(options.prefs.pp.profile, options.prefs.pp.csize);

    opts = cell(1, length(uipatdirs));
    for pat = 1:length(opts)

        % set patient specific options
        opts{pat} = options;
        opts{pat}.root = [fileparts(uipatdirs{pat}),filesep];
        [~, opts{pat}.patientname] = fileparts(uipatdirs{pat});
    end

    parfor pat = 1:length(uipatdirs)
        % run main function
        try
            ea_autocoord(opts{pat});
        catch
            warning([opts{pat}.patientname,' failed. Please run this patient again and adjust parameters. Moving on to next patient.']);
        end
    end

    delete(pp);
else
    switch cmd
        case 'export'
            % mark the lead task as exported job
            options.exportedJob = 1;
            ea_export(options);
        case 'run'
            for pat = 1:length(uipatdirs)

                % set patient specific options
                options.root = [fileparts(uipatdirs{pat}),filesep];
                [~, options.patientname] = fileparts(uipatdirs{pat});
                options.pat=pat;
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
    end
end

if needtoexit % close Matlab after job is done.
    exit
end
