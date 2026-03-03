function ea_connectomics_run_struct(app, event)
%EA_CONNECTOMICS_RUN_STRUCT
% Helper to run structural Lead-Connectome from the main lead_dbs app.
%
%   ea_connectomics_run_struct(app, event) is intended to be called from the
%   Connectomics tab "Run" button callback inside lead_dbs.mlapp
%   (App Designer). It:
%
%   - reads the current GUI state via ea_convertGUIDEArguments / handles
%   - uses the StructuralConnectivitySwitch in the Connectomics tab
%   - configures options.lc.* for STRUCTURAL ONLY
%   - calls ea_run('run', options) with options.leadprod = 'connectome'
%
% Usage inside lead_dbs.mlapp:
%
%   function RunButtonPushed(app, event)
%       ea_connectomics_run_struct(app, event);
%   end

    %#ok<*NASGU>  % hObject / eventdata not otherwise used here

    % Convert App Designer callback args to GUIDE-style
    [hObject, eventdata, handles] = ea_convertGUIDEArguments(app, event);

    % Show busy indicator in main GUI (if present)
    try
        ea_busyaction('on', handles.busystatus);
    catch
        % Fallback: ignore if busystatus handle not available
    end

    % Base options from main Lead-DBS GUI
    options = ea_handles2options(handles);
    options.uipatdirs = getappdata(handles.leadfigure, 'uipatdir');

    % Safety: require a selected patient
    if isempty(options.uipatdirs) || strcmp(options.uipatdirs{1}, 'No Patient Selected')
        try ea_busyaction('off', handles.busystatus); end
        ea_cprintf('CmdWinWarnings', 'No patient selected.\n');
        return;
    end

    % Derive root / patientname from first selected patient directory
    [options.root, options.patientname] = fileparts(options.uipatdirs{1});
    options.root = [options.root filesep];

    % ---------------------------
    %  Configure LC: STRUCTURAL
    % ---------------------------

    % Start from default LC prefs
    prefs = ea_prefs;

    lc.general.parcellation = prefs.lc.defaultParcellation;

    % Graph metrics defaults
    lc.graph.degree_centrality     = 0;
    lc.graph.eigenvector_centrality = 0;
    lc.graph.nodal_efficiency      = 0;
    lc.graph.struc_func_sim        = 0;
    lc.graph.fthresh               = nan;
    lc.graph.sthresh               = nan;

    % Functional completely OFF for now
    lc.func.compute_CM   = 0;
    lc.func.compute_GM   = 0;
    lc.func.prefs.TR     = 2;

    % Structural toggled by the upper switch in the Connectomics tab
    doStruc = false;
    try
        % App Designer switch has Values 'On' / 'Off'
        doStruc = strcmpi(app.StructuralConnectivitySwitch.Value, 'On');
    catch
        % If switch not found, default to ON for robustness
        doStruc = true;
    end

    lc.struc.ft.do        = doStruc;              % run fiber tracking
    lc.struc.ft.method    = 'ea_ft_gqi_yeh';      % default structural method
    lc.struc.ft.dsistudio.fiber_count = 200000;
    lc.struc.ft.normalize = doStruc;              % normalize fibers as well
    lc.struc.compute_CM   = doStruc;              % build DTI connectivity matrix
    lc.struc.compute_GM   = 0;                    % no graph metrics for now

    % NBS / advanced defaults
    lc.nbs.adv.method   = 1;
    lc.nbs.adv.compsize = 1;
    lc.nbs.adv.perm     = 5000;
    lc.nbs.adv.alpha    = 0.05;
    lc.nbs.adv.exch     = '';

    options.lc = lc;

    % Tell ea_autocoord to run the Lead-Connectome subroutine
    % (this triggers ea_perform_lc at the end of the pipeline)
    options.dolc = 1;

    % leadprod controls some other branches (e.g. mapper, predict).
    % We keep it as 'dbs' so that the standard DBS preprocessing
    % behaves as usual, but dolc=1 ensures the connectome part runs.
    options.leadprod = 'dbs';

    % Run structural Lead-Connectome pipeline for the current patient
    try
        ea_run('run', options);
    catch ME
        try ea_busyaction('off', handles.busystatus); end
        rethrow(ME);
    end

    try ea_busyaction('off', handles.busystatus); end
end

