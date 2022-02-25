function ea_genRawjson(subjDir)
% Generate rawimages.json based on the images in patient rawdata folder
% subjDir should be 'dataset/rawdata/sub-XXX', i.e., patient rawdata folder.

% Convert to cell
if ischar(subjDir)
    subjDir = {subjDir};
end

% Restore full path
subjDir = GetFullPath(subjDir);

% Remove filesep at the end
subjDir = regexprep(subjDir, ['\', filesep, '$'], '');

% Iterate subjDir
for i = 1:length(subjDir)
    rawdata = struct;

    niftiFiles = ea_regexpdir(subjDir{i}, '.*\.nii(\.gz)?$', 1, 'f');
    for f = 1:length(niftiFiles)
        parsed = parseBIDSFilePath(niftiFiles{f});
        session = parsed.ses; % preop, postop
        [~, type] = fileparts(parsed.dir); % anat, func, dwi
        if isfield(parsed, 'acq')
            suffix = [parsed.acq, '_', parsed.suffix]; % e.g., ax_T1w
        else
            suffix = parsed.suffix; % e.g., CT
        end
        [~, value] = ea_niifileparts(niftiFiles{f}); % File name without ext
        rawdata.(session).(type).(suffix) = value;
    end
    rawdata = orderfields(rawdata, {'preop', 'postop'}); % Re-order

    % Save rawimages.json
    jsonPath = regexprep(parsed.dir, ['rawdata\', filesep, '.*'], fullfile('derivatives', 'leaddbs', ['sub-', parsed.sub], 'prefs', ['sub-', parsed.sub, '_desc-rawimages.json']));
    fprintf('Generating rawimages.json:\n%s\n\n', jsonPath);
    savejson('', rawdata, jsonPath);
end
