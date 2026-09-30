function out = ea_runacpc(options)

    out = 0;

    preopFields = fieldnames(options.subj.preproc.anat.preop);

    preprocFiles = cellfun( ...
        @(x) options.subj.preproc.anat.preop.(x), ...
        preopFields, ...
        'UniformOutput', false);

    app = LeG_ACPCGUI('NiftiFile', preprocFiles{1});

    % Wait for the AC/PC GUI to explicitly finish.
    uiwait(app.MainFigure);

    % Record whether the image was actually saved/modified.
    if isvalid(app)
        out = app.SavedChanges;
        delete(app);
    end
end