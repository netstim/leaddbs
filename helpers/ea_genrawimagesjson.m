function rawImages = ea_genrawimagesjson(BIDSRoot, subjId)
% [Re-]generate rawimages json file in case it's not present in subject's
% derivatives folder

% Get all images
rawdataFolder = fullfile(GetFullPath(BIDSRoot), 'rawdata', ['sub-', subjId]);
preopNiftiFiles = ea_regexpdir([rawdataFolder, filesep, 'ses-preop'], '.*\.nii(\.gz)?$', 1, 'f');
postopNiftiFiles = ea_regexpdir([rawdataFolder, filesep, 'ses-postop'], '.*\.nii(\.gz)?$', 1, 'f');
niftiFiles = [preopNiftiFiles; postopNiftiFiles];

% Iterate all images
rawImages = struct;
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
    rawImages.(session).(type).(suffix) = value;
end

if isfield(rawImages, 'postop')
    rawImages = orderfields(rawImages, {'preop', 'postop'}); % Re-order
end

if isempty(fieldnames(rawImages))
    % Warn in case it's not a miniset
    if ~isfile(fullfile(GetFullPath(BIDSRoot), 'miniset.json'))
        ea_cprintf('CmdWinWarnings', 'No raw images found for "%s"!\n', subjId);
    end
    return;
end

ea_cprintf('CmdWinWarnings', 'Generating rawimages.json for "%s":\n', subjId);

% Get prefs folder
prefsFolder = fullfile(GetFullPath(BIDSRoot), 'derivatives', 'leaddbs', ['sub-', subjId], 'prefs');
ea_mkdir(prefsFolder);

% Save rawimages.json
jsonPath = fullfile(prefsFolder, ['sub-', subjId, '_desc-rawimages.json']);
ea_cprintf('CmdWinWarnings', '%s\n\n', jsonPath);
savejson('', rawImages, jsonPath);
