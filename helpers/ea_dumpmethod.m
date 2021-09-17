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
        else
            log.method.(modality) = options.coregmr.method;
        end
    case 'norm'
        log.method = options.normalize.method;
    case 'brainshift'
        log.method = options.scrf.mask;
end

% Save method log
savejson('', log, methodLogFile);
