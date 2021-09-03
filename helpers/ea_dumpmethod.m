function ea_dumpmethod(options, type)

methodLogFile = options.subj.(type).log.method;

% Load method log when it exists.
if isfile(methodLogFile)
    log = loadjson(methodLogFile);
end

% Set method
switch type
    case 'coreg'
        log.method.CT = options.coregct.method;
        log.method.MRI = options.coregmr.method;
    case 'norm'
        log.method = options.normalize.method;
end

% Save method log
savejson('', log, methodLogFile);
