function out = ea_runacpc(options)

    out = 0;

    preopFields = fieldnames(options.subj.preproc.anat.preop);

    preprocFiles = cellfun( ...
        @(x) options.subj.preproc.anat.preop.(x), ...
        preopFields, ...
        'UniformOutput', false);

    % ACPC log, following the existing Lead-DBS JSON pattern
    logDir = fullfile(options.subj.subjDir, 'preprocessing', 'log');
    logFile = fullfile(logDir, ...
    ['sub-', regexprep(options.subj.subjId, '^sub-', ''), '_desc-acpc.json']);

    ea_mkdir(logDir);

    % Preserve previous alignment status if a new run is canceled.
    json = struct('ACPCAligned', false);

    if isfile(logFile)
        previous = loadjson(logFile);
        if isfield(previous, 'ACPCAligned')
            json.ACPCAligned = logical(previous.ACPCAligned);
        end
    end

    app = LeG_ACPCGUI('NiftiFile', preprocFiles{1});

    % Wait for the ACPC GUI to finish.
    uiwait(app.MainFigure);

    if isvalid(app)
        out = app.SavedChanges;

        if out
            json.ACPCAligned = logical(app.ACPCAligned);
        end

        delete(app);
    end

    savejson('', json, logFile);

end
