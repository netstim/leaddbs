function ea_dumpmethod(options, type, modality)

if strcmp(type, 'coreg') && ~exist('modality', 'var')
    error('Missing modality when dumping coregistration method!');
end

methodLogFile = options.subj.(type).log.method;

ea_mkdir(fileparts(methodLogFile));

% Load method log when it exists.
if isfile(methodLogFile)
    log = loadjson(methodLogFile);
end

% Set method
switch type
    case 'coreg'
        if strcmp(modality, 'CT')
            log.method.CT = options.coregct.method;
            log.approval.CT = 0;
        else
            log.method.(modality) = options.coregmr.method;
            log.approval.(modality) = 0;
        end
    case 'norm'
        log.method = options.normalize.method;
        if isfield(log, 'settings') && ~isempty(fieldnames(log.settings))
            log.settings = struct;
        end
        if contains(log.method, 'Avants 2008')
            log.settings.preset = eval([options.prefs.machine.normsettings.ants_preset, '(''query'')']);
            log.settings.transform = options.prefs.machine.normsettings.ants_strategy;
            log.settings.metric = options.prefs.machine.normsettings.ants_metric;
            log.settings.subcorticalrefine = options.prefs.machine.normsettings.ants_scrf;
            log.settings.usefa = options.prefs.machine.normsettings.ants_usefa;
            log.settings.skullstripped = options.prefs.machine.normsettings.ants_skullstripped;
        elseif contains(log.method, 'Andersson 2010')
            log.settings.skullstrip = options.prefs.machine.normsettings.fsl_skullstrip;
        elseif contains(log.method, 'Schonecker 2009')
            if options.prefs.machine.normsettings.schoenecker_movim == 1
                log.settings.baseimage = 'pre-op';
            elseif options.prefs.machine.normsettings.schoenecker_movim == 2
                log.settings.baseimage = 'post-op';
            end
        elseif contains(log.method, 'Ashburner 2005')
            log.settings.regularization = options.prefs.machine.normsettings.spmnewseg_scalereg;
        else % Remove settings field when no settings available
            if isfield(log, 'settings')
                log = rmfield(log, 'settings');
            end
        end
        log.approval = 0;
    case 'brainshift'
        log.method = options.scrf.mask;
        log.approval = 0;
end

% Save method log
savejson('', log, methodLogFile);
