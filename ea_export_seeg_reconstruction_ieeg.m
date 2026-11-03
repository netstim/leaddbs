function [tsvFiles, coordsystemFiles] = ea_export_seeg_reconstruction_ieeg(reconstructionFile, opts)
%EA_EXPORT_SEEG_RECONSTRUCTION_IEEG Export native and MNI sEEG BIDS tables.
%   [TSV, JSON] = ea_export_seeg_reconstruction_ieeg(RECONSTRUCTIONFILE)
%   writes one *_electrodes.tsv and matching *_coordsystem.json per space:
%   reco.native.coords_mm and reco.mni.coords_mm. TSV and JSON are structs
%   with fields native and mni containing the resulting paths.
%
%   By default, the files go in an ieeg/ folder next to the reconstruction/
%   folder. For the example input at derivatives/leaddbs/sub-01/
%   reconstruction/sub-01_desc-reconstruction.mat, the destination is
%   derivatives/leaddbs/sub-01/ieeg/.
%
%   Name-value options:
%     TemplateSpace                Space of reco.mni (defaults to ea_getspace).
%     OutputDir                    Override destination ieeg/ folder.
%     Subject                      Override inferred sub- label.
%     Session                      Optional ses- label in output filenames.
%     NativeCoordinateSystem       BIDS keyword, default 'Other'.
%     NativeCoordinateDescription  Override native coordsystem description.
%     ContactAreaMM2               Scalar exposed contact area, default n/a.
%     Overwrite                    Replace existing exports, default false.
%
%   The MAT file does not identify the native-image origin/orientation or
%   contact surface area. Verify those details before using the export as
%   a complete iEEG localization record. For the MNI table, ensure the
%   selected TemplateSpace is the one used when the reconstruction was made.

arguments
    reconstructionFile (1,1) string {mustBeFile}
    opts.TemplateSpace (1,1) string = ""
    opts.OutputDir (1,1) string = ""
    opts.Subject (1,1) string = ""
    opts.Session (1,1) string = ""
    opts.NativeCoordinateSystem (1,1) string {mustBeMember(opts.NativeCoordinateSystem, ["Other", "ScanRAS", "ACPC"])} = "Other"
    opts.NativeCoordinateDescription (1,1) string = ""
    opts.ContactAreaMM2 (1,1) double = NaN
    opts.Overwrite (1,1) logical = false
end

if ~isnan(opts.ContactAreaMM2) && (~isfinite(opts.ContactAreaMM2) || opts.ContactAreaMM2 <= 0)
    error('sEEG:BIDS:InvalidSize', 'ContactAreaMM2 must be positive or NaN.');
end

data = load(reconstructionFile, 'reco');
if ~isfield(data, 'reco') || ~isstruct(data.reco) || ~isscalar(data.reco)
    error('sEEG:BIDS:InvalidReconstruction', 'The MAT file must contain a scalar reco struct.');
end
reco = data.reco;
if ~isfield(reco, 'props') || ~isstruct(reco.props) || isempty(reco.props)
    error('sEEG:BIDS:MissingProperties', 'reco.props is missing or empty.');
end
spaces = {'native', 'mni'};
coordinates = struct;
for spaceIndex = 1:numel(spaces)
    space = spaces{spaceIndex};
    if ~isfield(reco, space) || ~isstruct(reco.(space)) || ...
            ~isfield(reco.(space), 'coords_mm')
        error('sEEG:BIDS:MissingCoordinates', 'reco.%s.coords_mm is missing.', space);
    end
    coordinateCells = reco.(space).coords_mm;
    if isnumeric(coordinateCells) && isscalar(reco.props)
        coordinateCells = {coordinateCells};
    end
    if ~iscell(coordinateCells) || numel(coordinateCells) ~= numel(reco.props)
        error('sEEG:BIDS:InvalidCoordinates', ...
            'reco.%s.coords_mm must have one array per electrode.', space);
    end
    coordinates.(space) = coordinateCells;
end

[matDir, matBase] = fileparts(char(reconstructionFile));
subject = char(opts.Subject);
if isempty(subject)
    subjectMatch = regexp(matBase, '^sub-([A-Za-z0-9]+)(?:_.*)?$', 'tokens', 'once');
    if isempty(subjectMatch)
        [subjectDir, folderName] = fileparts(matDir);
        if strcmp(folderName, 'reconstruction')
            [~, folderName] = fileparts(subjectDir);
        end
        subjectMatch = regexp(folderName, '^sub-([A-Za-z0-9]+)$', 'tokens', 'once');
    end
    if isempty(subjectMatch)
        error('sEEG:BIDS:MissingSubject', ...
            'Cannot infer the BIDS subject label; pass Subject explicitly.');
    end
    subject = subjectMatch{1};
end
subject = regexprep(subject, '^sub-', '');
assertBidsLabel(subject, 'Subject');

if strlength(opts.OutputDir) == 0
    [parentDir, folderName] = fileparts(matDir);
    if strcmp(folderName, 'reconstruction')
        outputDir = fullfile(parentDir, 'ieeg');
    else
        outputDir = fullfile(matDir, 'ieeg');
    end
else
    outputDir = char(opts.OutputDir);
end
[sessionParent, outputFolderName] = fileparts(outputDir);
if ~strcmp(outputFolderName, 'ieeg')
    error('sEEG:BIDS:InvalidOutputDir', 'OutputDir must be an ieeg/ folder.');
end
session = regexprep(char(opts.Session), '^ses-', '');
if isempty(session)
    [~, parentName] = fileparts(sessionParent);
    if startsWith(parentName, 'ses-')
        session = char(extractAfter(parentName, 4));
    end
end
if ~isempty(session)
    assertBidsLabel(session, 'Session');
end

templateSpace = char(opts.TemplateSpace);
if isempty(templateSpace)
    if exist('ea_getspace', 'file') ~= 2
        error('sEEG:BIDS:MissingTemplateSpace', ...
            'Pass TemplateSpace when ea_getspace is unavailable.');
    end
    templateSpace = char(ea_getspace);
    warning('sEEG:BIDS:VerifyTemplateSpace', ...
        'Using current Lead-DBS space %s; confirm it matches this reconstruction.', templateSpace);
end
templateLabel = regexprep(templateSpace, '[^A-Za-z0-9]', '');
assertBidsLabel(templateLabel, 'TemplateSpace');

% Determine native coordinate system from the ACPC log
logFile = fullfile(sessionParent, 'preprocessing', 'log', ...
    ['sub-', subject, '_desc-acpc.json']);

% Default to ScanRAS when ACPC has not been performed
nativeSystem = 'ScanRAS';

if isfile(logFile)
    acpcLog = loadjson(logFile);

    if isfield(acpcLog, 'ACPCAligned') && acpcLog.ACPCAligned
        nativeSystem = 'ACPC';
    end
end

baseName = ['sub-', subject];
if ~isempty(session)
    baseName = [baseName, '_ses-', session];
end

tsvFiles = struct( ...
    'native', fullfile(outputDir, [baseName, '_space-', nativeSystem, '_electrodes.tsv']), ...
    'mni', fullfile(outputDir, [baseName, '_space-', templateLabel, '_electrodes.tsv']));

coordsystemFiles = struct( ...
    'native', fullfile(outputDir, [baseName, '_space-', nativeSystem, '_coordsystem.json']), ...
    'mni', fullfile(outputDir, [baseName, '_space-', templateLabel, '_coordsystem.json']));

destinations = {tsvFiles.native, coordsystemFiles.native, ...
    tsvFiles.mni, coordsystemFiles.mni};
if ~opts.Overwrite && any(cellfun(@isfile, destinations))
    error('sEEG:BIDS:Exists', ...
        'An iEEG export already exists. Pass Overwrite=true to replace it.');
end

[nativeText, nativeNames] = buildTsv(reco.props, coordinates.native, opts.ContactAreaMM2);
[mniText, mniNames] = buildTsv(reco.props, coordinates.mni, opts.ContactAreaMM2);
if ~isequal(nativeNames, mniNames)
    error('sEEG:BIDS:ContactMismatch', 'Native and MNI contact names do not match.');
end


nativeDescription = char(opts.NativeCoordinateDescription);
if isempty(nativeDescription)
    nativeDescription = ['Lead-DBS native reconstruction coordinates in millimeters. ', ...
        'Confirm origin and orientation against the subject-specific anatomical ', ...
        'image used for reconstruction.'];
end
if isempty(nativeSystem)
    error('sEEG:BIDS:InvalidCoordinateSystem', 'NativeCoordinateSystem cannot be empty.');
end
nativeMetadata = makeCoordinateMetadata(nativeSystem, nativeDescription, 'native');

standardSpaces = {'MNI152NLin2009aSym', 'MNI152NLin2009bSym', ...
    'MNI152NLin2009cSym', 'MNI152NLin2009aAsym', ...
    'MNI152NLin2009bAsym', 'MNI152NLin2009cAsym', ...
    'MNI152NLin6Sym', 'MNI152NLin6Asym', 'MNI152Lin', ...
    'MNI305', 'MNIColin27', 'OASIS30Atropos'};
if ismember(templateSpace, standardSpaces)
    mniSystem = templateSpace;
else
    mniSystem = 'Other';
end
mniDescription = sprintf('Lead-DBS template-space reconstruction coordinates in %s.', templateSpace);
mniMetadata = makeCoordinateMetadata(mniSystem, mniDescription, 'mni');

if ~isfolder(outputDir)
    [made, message] = mkdir(outputDir);
    if ~made
        error('sEEG:BIDS:CreateDirectoryFailed', '%s', message);
    end
end
writeTextAtomically(tsvFiles.native, nativeText);
writeTextAtomically(coordsystemFiles.native, [jsonencode(nativeMetadata), newline]);
writeTextAtomically(tsvFiles.mni, mniText);
writeTextAtomically(coordsystemFiles.mni, [jsonencode(mniMetadata), newline]);
fprintf('Exported %d contacts in native and MNI spaces to %s\n', numel(nativeNames), outputDir);
end

function [text, names] = buildTsv(props, coordinateCells, contactArea)
lines = {sprintf('name\tx\ty\tz\tsize\tgroup\ttype')};
names = strings(0,1);
for electrodeIndex = 1:numel(props)
    prop = props(electrodeIndex);
    xyzRows = coordinateCells{electrodeIndex};
    if isempty(xyzRows)
        xyzRows = zeros(0,3);
    end
    if ~isnumeric(xyzRows) || ~isreal(xyzRows) || ~ismatrix(xyzRows) || size(xyzRows,2) ~= 3
        error('sEEG:BIDS:InvalidCoordinates', ...
            'Electrode %d must have an N-by-3 numeric coordinate array.', electrodeIndex);
    end

    group = sprintf('lead%02d', electrodeIndex);
    if isfield(prop, 'elname') && ~isempty(prop.elname)
        group = char(string(prop.elname));
    end
    assertTsvText(group, sprintf('electrode %d group', electrodeIndex));

    labels = {};
    if isfield(prop, 'labels') && ~isempty(prop.labels)
        labels = prop.labels;
        if iscell(labels) && isscalar(labels) && iscell(labels{1})
            labels = labels{1};
        end
        if isstring(labels) || ischar(labels)
            labels = cellstr(labels);
        end
        if ~iscell(labels)
            error('sEEG:BIDS:InvalidLabels', 'Electrode %d labels are not text.', electrodeIndex);
        end
        labels = labels(:);
    end
    if isempty(labels)
        labels = arrayfun(@(i) sprintf('%s%d', group, i), ...
            1:size(xyzRows,1), 'UniformOutput', false)';
    end
    if ~isempty(xyzRows) && numel(labels) ~= size(xyzRows,1)
        error('sEEG:BIDS:LabelMismatch', ...
            'Electrode %d has %d labels but %d coordinate rows.', ...
            electrodeIndex, numel(labels), size(xyzRows,1));
    end

    for contactIndex = 1:numel(labels)
        name = char(string(labels{contactIndex}));
        assertTsvText(name, sprintf('electrode %d contact %d name', electrodeIndex, contactIndex));
        if isempty(name)
            error('sEEG:BIDS:InvalidLabels', 'Contact names cannot be empty.');
        end
        names(end+1,1) = string(name); %#ok<AGROW>
        xyz = {'n/a', 'n/a', 'n/a'};
        if size(xyzRows,1) >= contactIndex && all(isfinite(xyzRows(contactIndex,:)))
            xyz = arrayfun(@(v) sprintf('%.10g', v), xyzRows(contactIndex,:), ...
                'UniformOutput', false);
        end
        area = 'n/a';
        if isfinite(contactArea)
            area = sprintf('%.10g', contactArea);
        end
        lines{end+1,1} = strjoin([{name}, xyz, {area, group, 'depth'}], sprintf('\t')); %#ok<AGROW>
    end
end
if isempty(names)
    error('sEEG:BIDS:NoContacts', 'The reconstruction contains no contacts.');
end
if numel(unique(names)) ~= numel(names)
    error('sEEG:BIDS:DuplicateNames', 'Contact names must be unique across electrodes.');
end
text = [strjoin(lines, newline), newline];
end

function metadata = makeCoordinateMetadata(coordinateSystem, description, sourceSpace)
metadata = struct( ...
    'iEEGCoordinateSystem', coordinateSystem, ...
    'iEEGCoordinateUnits', 'mm', ...
    'iEEGCoordinateSystemDescription', description, ...
    'iEEGCoordinateProcessingDescription', ...
        sprintf('Contact centers exported from Lead-DBS reco.%s.coords_mm.', sourceSpace));
end

function assertBidsLabel(value, name)
if isempty(regexp(value, '^[A-Za-z0-9]+$', 'once'))
    error('sEEG:BIDS:InvalidLabel', '%s must contain only letters and digits.', name);
end
end

function assertTsvText(value, name)
if any(ismember(value, [char(9), newline, char(13)]))
    error('sEEG:BIDS:InvalidText', '%s contains a tab or newline.', name);
end
end

function writeTextAtomically(destination, contents)
temporaryFile = [tempname(fileparts(destination)), '.tmp'];
[fileID, message] = fopen(temporaryFile, 'w', 'n', 'UTF-8');
if fileID < 0
    error('sEEG:BIDS:WriteFailed', 'Cannot open temporary file: %s', message);
end
fileCleanup = onCleanup(@() fclose(fileID));
fprintf(fileID, '%s', contents);
clear fileCleanup
[moved, message] = movefile(temporaryFile, destination, 'f');
if ~moved
    if isfile(temporaryFile)
        delete(temporaryFile);
    end
    error('sEEG:BIDS:WriteFailed', 'Cannot write %s: %s', destination, message);
end
end
