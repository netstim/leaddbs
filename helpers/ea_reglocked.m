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
    case 'brainshift'
        key = 'brainshift';
    case 'normalization'
        key = 'norm';
end

% Check log file
logFile = options.subj.(key).log.method;
if ~isfile(logFile) % log file doesn't exist
    return;
else
    % Load method log
    json = loadjson(logFile);

    % Extract image modality
    modality = regexp(imagePath, '(?<=_)([a-zA-Z0-9]+)(?=\.nii(\.gz)?$)', 'match', 'once');

    % Check approval status
    if ~isfield(json, 'approval') || ~isfield(json.approval.(modality))
        return;
    else
        locked = json.approval.(modality);
    end
end
