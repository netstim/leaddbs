function locked = ea_reglocked(options, imagePath)
% Check if registration is locked (has been approved).
% Return 1 if so, otherwise return 0.

locked = 0;

if isfield(options, 'overwriteapproved') && options.overwriteapproved
    return
end

% Check pipeline keyword
[~, pipeline] = fileparts(fileparts(fileparts(imagePath)));
switch pipeline
    case 'coregistration'
        key = 'coreg';
    case 'normalization'
        key = 'norm';
    case 'brainshift'
        key = 'brainshift';
end

% Check log file
logFile = options.subj.(key).log.method;
if isfile(logFile)
    % Load method log
    json = loadjson(logFile);

    % Extract image modality
    modality = ea_getmodality(imagePath);

    % Check approval status
    if eval('isfield(json.approval, modality)', '0')
        locked = json.approval.(modality);
    end
end
