function out = ea_runacpc(options)
    preproc_files = cellfun(@(x) options.subj.preproc.anat.preop.(x), fieldnames(options.subj.preproc.anat.preop), 'uni', 0);
    app = LeG_ACPCGUI('NiftiFile', preproc_files{1});   % and optionally:
    out = 1;
end