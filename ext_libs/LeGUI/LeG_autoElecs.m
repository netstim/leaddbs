function [WC, T, debug] = LeG_autoElecs(app, varargin)
%LEG_AUTOELECS Detect intracranial electrode contacts from a postoperative CT.
%
% This is the MATLAB-only runtime for the frozen population contact model.
% It does not call Python and does not use reference contacts, electrode
% counts, electrode labels, or patient identity during detection.
%
% Minimal Lead-DBS invocation:
%   p = struct('subject', 'sub-0d5e8e', ...
%       'leadRoot', '/path/to/derivatives/leaddbs');
%   [WC, T, debug] = LeG_autoElecs(app, p);
%
% Direct CT invocation:
%   [WC, T, debug] = LeG_autoElecs(app, 'ctPath', '/path/postop_CT.nii');
%
% Principal operating parameters:
%   detectionThreshold  classifier probability threshold (default 0.55)
%   pointSeparationMm   minimum distance between outputs (default 2.5 mm)
%   brainMarginMm       tissue-mask expansion (default 30 mm)
%   intensityPercentile bright candidate threshold (default 99.0)
%   sourcePercentile    very-bright indicator threshold (default 99.8)
%   gridSpacingMm       candidate-grid spacing (default 1.0 mm)
%   maxCandidates       maximum candidates scored (default 60000)

started = tic;
options = localParseOptions(app, varargin{:});
model = localLoadModel(options.modelPath);
options = localApplyModelDefaults(options, model);
localValidateDependencies();

subject = localSubject(options.subject, app);
[ctPath, anatFolder] = localResolveCTPath(options.ctPath, options.leadRoot, subject, app);
localPrint(options, '[AutoElec] MATLAB population detector: %s\n', subject);
localPrint(options, '[AutoElec] CT: %s\n', ctPath);

stage = tic;
[rawCT, ctAffine] = localReadNifti(ctPath);
finiteCT = isfinite(rawCT);
if ~any(finiteCT(:))
    error('LeG_autoElecs:EmptyCT', 'The CT contains no finite voxels.');
end
[tissue, maskPaths] = localLoadTissueMasks( ...
    anatFolder, options.maskPaths, size(rawCT), ctAffine, options);
voxelSizes = sqrt(sum(ctAffine(1:3, 1:3).^2, 1));
[support, headCenter, brainMaskAvailable] = localBrainSupport( ...
    tissue, finiteCT, ctAffine, voxelSizes, options);
loadSeconds = toc(stage);

stage = tic;
[pool, candidateMeta] = localBuildCandidates( ...
    rawCT, finiteCT, support, tissue, ctAffine, options);
candidateSeconds = toc(stage);
localPrint(options, ...
    '[AutoElec] Candidates: %d at CT p%.3g (threshold %.3f)\n', ...
    size(pool.points, 1), options.intensityPercentile, candidateMeta.intensityThreshold);
if isempty(pool.points)
    WC = zeros(0, 3);
    T = NaN;
    debug = localEmptyDebug(subject, ctPath, maskPaths, options, model, started);
    return;
end

stage = tic;
[topology, pool] = localPopulationFeatures(pool, options.neighborBatchSize);
normalizedCT = localNormalizeCT(rawCT, finiteCT);
appearance = localDenseFeatures(pool, normalizedCT, ctAffine, headCenter);
features = [topology, appearance];
featureSeconds = toc(stage);
featureNames = localFeatureNames();
if size(features, 2) ~= double(model.featureCount)
    error('LeG_autoElecs:FeatureCount', ...
        'Generated %d features, but the model requires %d.', ...
        size(features, 2), double(model.featureCount));
end
if ~isequal(featureNames(:), localCellstr(model.featureNames(:)))
    error('LeG_autoElecs:FeatureSchema', ...
        'The generated feature order does not match the model artifact.');
end

stage = tic;
scores = localPredictLightGBM(features, model);
eligible = find(scores >= options.detectionThreshold);
selectedLocal = localSpatialNms( ...
    pool.points(eligible, :), scores(eligible), options.pointSeparationMm);
selected = eligible(selectedLocal);
predictionMm = pool.points(selected, :);
selectedScores = scores(selected);
inferenceSeconds = toc(stage);

[WC, outputTransform] = localOutputVoxels(app, predictionMm, ctAffine);
T = NaN;
debug = struct;
debug.mode = 'population_contact_classifier_matlab';
debug.subject = subject;
debug.final_mm = predictionMm;
debug.wc = WC;
debug.selected_scores = selectedScores;
debug.selected_candidate_indices = selected;
debug.candidate_mm = pool.points;
debug.candidate_scores = scores;
debug.reference_loaded = false;
debug.contact_count_used = false;
debug.patient_identity_used = false;
debug.ct_path = ctPath;
debug.mask_paths = maskPaths;
debug.brain_mask_available = brainMaskAvailable;
debug.model_path = options.modelPath;
debug.model_training_scope = localText(model.trainingScope);
debug.output_transform = outputTransform;
debug.parameters = options;
debug.candidate_metadata = candidateMeta;
debug.feature_names = featureNames;
if options.returnFeatures
    debug.features = features;
end
debug.timings = struct('load_seconds', loadSeconds, ...
    'candidate_seconds', candidateSeconds, ...
    'feature_seconds', featureSeconds, ...
    'inference_seconds', inferenceSeconds, ...
    'total_seconds', toc(started));

if ~isempty(options.savePath)
    localSaveResult(options.savePath, predictionMm, selectedScores, debug);
end
localPrint(options, ...
    '[AutoElec] Returned %d contacts from %d candidates in %.1f s.\n', ...
    size(WC, 1), size(pool.points, 1), toc(started));
end

function options = localParseOptions(app, varargin)
here = fileparts(mfilename('fullpath'));
options = struct;
options.subject = '';
options.leadRoot = '';
options.ctPath = '';
options.modelPath = fullfile(here, 'LeG_autoElecs_population_model.mat');
options.maskPaths = {};
options.detectionThreshold = [];
options.pointSeparationMm = [];
options.brainMarginMm = [];
options.intensityPercentile = [];
options.sourcePercentile = [];
options.gridSpacingMm = [];
options.maxCandidates = [];
options.useBrainMask = true;
options.returnFeatures = false;
options.verbose = true;
options.savePath = '';
options.neighborBatchSize = 2000;

detect = localReadApp(app, 'detect', struct);
if isstruct(detect)
    options = localMergeRecognized(options, detect);
end
if isempty(varargin)
    return;
end
if isscalar(varargin) && isstruct(varargin{1})
    options = localMergeRecognized(options, varargin{1});
elseif mod(numel(varargin), 2) == 0
    supplied = struct;
    for index = 1:2:numel(varargin)
        name = varargin{index};
        if ~(ischar(name) || isstring(name))
            error('LeG_autoElecs:InvalidOption', 'Option names must be text.');
        end
        supplied.(char(name)) = varargin{index + 1};
    end
    options = localMergeRecognized(options, supplied);
else
    error('LeG_autoElecs:InvalidOptions', ...
        'Pass one settings struct or name-value pairs.');
end
end

function options = localMergeRecognized(options, supplied)
fields = fieldnames(supplied);
for index = 1:numel(fields)
    value = supplied.(fields{index});
    key = lower(regexprep(fields{index}, '[^a-zA-Z0-9]', ''));
    switch key
        case {'subject', 'subjectid'}
            options.subject = localText(value);
        case {'leadroot', 'leaddbsroot'}
            options.leadRoot = localText(value);
        case {'ctpath', 'postopct'}
            options.ctPath = localText(value);
        case {'modelpath', 'populationmodelpath', 'pointmodelpath'}
            options.modelPath = localText(value);
        case {'maskpaths', 'tissuemaskpaths'}
            options.maskPaths = cellstr(string(value));
        case {'detectionthreshold', 'threshold'}
            options.detectionThreshold = double(value(1));
        case {'pointseparationmm', 'separationmm'}
            options.pointSeparationMm = double(value(1));
        case 'brainmarginmm'
            options.brainMarginMm = double(value(1));
        case {'intensitypercentile', 'candidatepercentile'}
            options.intensityPercentile = double(value(1));
        case 'sourcepercentile'
            options.sourcePercentile = double(value(1));
        case 'gridspacingmm'
            options.gridSpacingMm = double(value(1));
        case 'maxcandidates'
            options.maxCandidates = round(double(value(1)));
        case 'usebrainmask'
            options.useBrainMask = logical(value(1));
        case 'returnfeatures'
            options.returnFeatures = logical(value(1));
        case 'verbose'
            options.verbose = logical(value(1));
        case {'savepath', 'outputpath'}
            options.savePath = localText(value);
        case 'neighborbatchsize'
            options.neighborBatchSize = round(double(value(1)));
        otherwise
            % Legacy GUI and Python-wrapper settings are intentionally ignored.
    end
end
end

function model = localLoadModel(path)
if ~isfile(path)
    error('LeG_autoElecs:MissingModel', 'Model file not found: %s', path);
end
loaded = load(path, 'model');
if ~isfield(loaded, 'model') || ~isstruct(loaded.model)
    error('LeG_autoElecs:InvalidModel', '%s does not contain a model struct.', path);
end
model = loaded.model;
required = {'featureNames', 'featureCount', 'roots', 'featureIndex', ...
    'threshold', 'leftChild', 'rightChild', 'defaultLeft', 'isLeaf', ...
    'leafValue', 'decisionThreshold', 'separationMm'};
for index = 1:numel(required)
    if ~isfield(model, required{index})
        error('LeG_autoElecs:InvalidModel', ...
            'Model field %s is missing from %s.', required{index}, path);
    end
end
end

function options = localApplyModelDefaults(options, model)
defaults = {
    'detectionThreshold', 'decisionThreshold';
    'pointSeparationMm', 'separationMm';
    'brainMarginMm', 'brainMarginMm';
    'intensityPercentile', 'intensityPercentile';
    'sourcePercentile', 'sourcePercentile';
    'gridSpacingMm', 'gridSpacingMm';
    'maxCandidates', 'maxCandidates'};
for index = 1:size(defaults, 1)
    optionName = defaults{index, 1};
    modelName = defaults{index, 2};
    if isempty(options.(optionName))
        options.(optionName) = double(model.(modelName));
    end
end
if ~isfinite(options.detectionThreshold) || options.detectionThreshold < 0 || ...
        options.detectionThreshold > 1
    error('LeG_autoElecs:Threshold', 'detectionThreshold must be in [0, 1].');
end
positive = {'pointSeparationMm', 'brainMarginMm', 'gridSpacingMm', ...
    'maxCandidates', 'neighborBatchSize'};
for index = 1:numel(positive)
    if ~isfinite(options.(positive{index})) || options.(positive{index}) <= 0
        error('LeG_autoElecs:OptionRange', '%s must be positive.', positive{index});
    end
end
end

function localValidateDependencies()
required = {'niftiinfo', 'niftiread', 'knnsearch', 'rangesearch', ...
    'KDTreeSearcher', 'bwdist', 'imresize3', 'prctile'};
missing = required(cellfun(@(name) exist(name, 'file') == 0 && ...
    exist(name, 'class') == 0, required));
if ~isempty(missing)
    error('LeG_autoElecs:MissingToolbox', ...
        'Required MATLAB functions are unavailable: %s.', strjoin(missing, ', '));
end
end

function subject = localSubject(requested, app)
subject = strtrim(localText(requested));
if isempty(subject)
    names = {'PatientIDStr', 'subject', 'Subject', 'PatientID'};
    for index = 1:numel(names)
        subject = strtrim(localText(localReadApp(app, names{index}, '')));
        if ~isempty(subject)
            break;
        end
    end
end
if isempty(subject)
    subject = 'sub-unknown';
elseif ~startsWith(subject, 'sub-')
    subject = ['sub-' subject];
end
end

function [ctPath, anatFolder] = localResolveCTPath(requested, leadRoot, subject, app)
ctPath = strtrim(localText(requested));
if isempty(ctPath)
    try
        if ~isempty(app.CTInfo) && isfield(app.CTInfo, 'fname')
            ctPath = app.CTInfo.fname;
        end
    catch
    end
end

if isempty(ctPath)
    appNames = {'CTPath', 'ctPath', 'PostopCTPath', 'CTFile'};
    for index = 1:numel(appNames)
        ctPath = strtrim(localText(localReadApp(app, appNames{index}, '')));
        if ~isempty(ctPath)
            break;
        end
    end
end
if isempty(ctPath) && ~isempty(leadRoot) && ~strcmp(subject, 'sub-unknown')
    roots = {leadRoot, fullfile(leadRoot, 'derivatives', 'leaddbs')};
    for rootIndex = 1:numel(roots)
        folder = fullfile(roots{rootIndex}, subject, 'coregistration', 'anat');
        hits = [dir(fullfile(folder, ...
            '*_ses-postop_space-anchorNative_desc-preproc_CT.nii')); ...
            dir(fullfile(folder, ...
            '*_ses-postop_space-anchorNative_desc-preproc_CT.nii.gz'))];
        hits = hits(~startsWith({hits.name}, '._'));
        if ~isempty(hits)
            [~, order] = sort({hits.name});
            hits = hits(order);
            ctPath = fullfile(hits(1).folder, hits(1).name);
            break;
        end
    end
end
if isempty(ctPath)
    error('LeG_autoElecs:MissingCT', ...
        'Provide ctPath, or provide subject and leadRoot.');
end
if ~isfile(ctPath)
    error('LeG_autoElecs:MissingCT', 'CT file not found: %s', ctPath);
end
anatFolder = fileparts(ctPath);
end

function [volume, affine, info] = localReadNifti(path)
info = niftiinfo(path);
stored = niftiread(info);
slope = double(info.MultiplicativeScaling);
offset = double(info.AdditiveOffset);
if isempty(slope) || ~isfinite(slope) || slope == 0
    slope = 1;
end
if isempty(offset) || ~isfinite(offset)
    offset = 0;
end
volume = single(stored) .* single(slope) + single(offset);
affine = double(info.Transform.T');
if ~isequal(size(affine), [4, 4]) || abs(affine(4, 4)) < eps
    error('LeG_autoElecs:InvalidAffine', 'Invalid NIfTI transform in %s.', path);
end
end

function [tissue, paths] = localLoadTissueMasks( ...
        anatFolder, requested, ctSize, ctAffine, options)
paths = requested;
if isempty(paths)
    paths = {};
    labels = {'GM', 'WM', 'CSF'};
    for labelIndex = 1:numel(labels)
        hits = [dir(fullfile(anatFolder, ...
            ['*label-' labels{labelIndex} '*_mask.nii'])); ...
            dir(fullfile(anatFolder, ...
            ['*label-' labels{labelIndex} '*_mask.nii.gz']))];
        hits = hits(~startsWith({hits.name}, '._'));
        for hitIndex = 1:numel(hits)
            paths{end + 1, 1} = fullfile(hits(hitIndex).folder, hits(hitIndex).name); %#ok<AGROW>
        end
    end
end
paths = unique(cellstr(string(paths)), 'stable');
tissue = zeros(ctSize, 'single');
validPaths = {};
for index = 1:numel(paths)
    path = paths{index};
    if ~isfile(path)
        warning('LeG_autoElecs:MissingMask', 'Ignoring missing tissue mask: %s', path);
        continue;
    end
    [mask, maskAffine] = localReadNifti(path);
    if ~isequal(size(mask), ctSize) || max(abs(maskAffine(:) - ctAffine(:))) > 1e-3
        localPrint(options, '[AutoElec] Resampling tissue mask: %s\n', path);
        mask = localResampleVolume(mask, maskAffine, ctSize, ctAffine);
    end
    tissue = max(tissue, mask);
    validPaths{end + 1, 1} = path; %#ok<AGROW>
end
paths = validPaths;
end

function output = localResampleVolume(source, sourceAffine, targetSize, targetAffine)
output = zeros(targetSize, 'single');
blockDepth = 8;
for zStart = 1:blockDepth:targetSize(3)
    zEnd = min(zStart + blockDepth - 1, targetSize(3));
    [x, y, z] = ndgrid(0:targetSize(1) - 1, ...
        0:targetSize(2) - 1, zStart - 1:zEnd - 1);
    voxels = [x(:), y(:), z(:)];
    world = localApplyAffine(targetAffine, voxels);
    values = localSampleWorld(source, sourceAffine, world);
    output(:, :, zStart:zEnd) = reshape(values, ...
        [targetSize(1), targetSize(2), zEnd - zStart + 1]);
end
end

function [support, centerMm, available] = localBrainSupport( ...
        tissue, finiteCT, affine, voxelSizes, options)
brain = tissue > 0.02;
available = nnz(brain) >= 100;
if ~available
    support = finiteCT;
    centerVoxel = (double(size(tissue)) - 1) ./ 2;
    centerMm = localApplyAffine(affine, centerVoxel);
    warning('LeG_autoElecs:NoTissueMasks', ...
        ['No usable Lead-DBS GM/WM/CSF masks were found. Detection will run ' ...
        'without the trained skull-suppression support mask.']);
    return;
end
centerVoxel = localMaskCenter(brain);
centerMm = localApplyAffine(affine, centerVoxel);
if ~options.useBrainMask
    support = finiteCT;
    return;
end

anisotropy = max(voxelSizes) / max(min(voxelSizes), eps);
if anisotropy <= 1.1
    distanceMm = bwdist(brain) .* mean(voxelSizes);
    support = distanceMm <= options.brainMarginMm;
else
    targetSpacing = 1.0;
    isotropicSize = max(round(double(size(brain)) .* voxelSizes ./ targetSpacing), 1);
    isotropicBrain = imresize3(brain, isotropicSize, 'nearest');
    isotropicSupport = bwdist(isotropicBrain) .* targetSpacing <= options.brainMarginMm;
    support = imresize3(isotropicSupport, size(brain), 'nearest');
end
support = logical(support) & finiteCT;
end

function center = localMaskCenter(mask)
count = double(nnz(mask));
profileX = double(squeeze(sum(sum(mask, 2), 3)));
profileY = double(squeeze(sum(sum(mask, 1), 3)));
profileZ = double(squeeze(sum(sum(mask, 1), 2)));
center = [dot((0:numel(profileX) - 1)', profileX), ...
    dot((0:numel(profileY) - 1)', profileY), ...
    dot((0:numel(profileZ) - 1)', profileZ)] ./ count;
end

function [pool, metadata] = localBuildCandidates( ...
        raw, finite, support, tissue, affine, options)
population = raw(finite & support);
if isempty(population)
    error('LeG_autoElecs:EmptySupport', 'The CT support mask contains no finite voxels.');
end
percentiles = prctile(population, ...
    [50, options.intensityPercentile, options.sourcePercentile, 99.9]);
p50 = double(percentiles(1));
intensityThreshold = double(percentiles(2));
sourceThreshold = double(percentiles(3));
p999 = double(percentiles(4));
bright = finite & support & raw >= intensityThreshold;
linear = find(bright);
if isempty(linear)
    pool = localEmptyPool;
    metadata = struct('intensityThreshold', intensityThreshold, ...
        'sourceIntensityThreshold', sourceThreshold, ...
        'supportVoxels', nnz(support), 'brightVoxels', 0);
    return;
end
[x, y, z] = ind2sub(size(raw), linear);
voxels = double([x - 1, y - 1, z - 1]);
intensity = double(raw(linear));
points = localApplyAffine(affine, voxels);
bins = int32(floor(points ./ options.gridSpacingMm));
% The training candidate builder sorted ascending and then reversed the
% index vector. Reproduce that tie direction explicitly.
[~, order] = sort(intensity, 'ascend');
order = flipud(order);
[~, first] = unique(bins(order, :), 'rows', 'stable');
keep = order(first);
points = points(keep, :);
voxels = voxels(keep, :);
intensity = intensity(keep);
if size(points, 1) > options.maxCandidates
    [~, strongest] = sort(intensity, 'ascend');
    strongest = flipud(strongest);
    strongest = strongest(1:options.maxCandidates);
    points = points(strongest, :);
    voxels = voxels(strongest, :);
    intensity = intensity(strongest);
end
tissueAtPoint = double(localSampleWorld(tissue, affine, points));
source = double(intensity >= sourceThreshold);
base = min(max((intensity - p50) ./ max(p999 - p50, 1e-6), 0), 4);
pool = struct('points', points, 'voxels', voxels, ...
    'baseScore', base(:), 'sourceBase', source(:), ...
    'tissue', tissueAtPoint(:), 'nn1', [], 'count5', []);
metadata = struct('intensityThreshold', intensityThreshold, ...
    'sourceIntensityThreshold', sourceThreshold, 'p50', p50, 'p999', p999, ...
    'supportVoxels', nnz(support), 'brightVoxels', numel(linear), ...
    'candidateCount', size(points, 1), 'samplingMode', 'bright_grid');
end

function pool = localEmptyPool
pool = struct('points', zeros(0, 3), 'voxels', zeros(0, 3), ...
    'baseScore', zeros(0, 1), 'sourceBase', zeros(0, 1), ...
    'tissue', zeros(0, 1), 'nn1', zeros(0, 1), 'count5', zeros(0, 1));
end

function normalized = localNormalizeCT(raw, finite)
percentiles = prctile(raw(finite), [95, 99.9]);
denominator = max(double(percentiles(2) - percentiles(1)), 1e-6);
normalized = max((raw - percentiles(1)) ./ single(denominator), 0);
normalized(~finite) = 0;
end

function [features, pool] = localPopulationFeatures(pool, batchSize)
points = pool.points;
n = size(points, 1);
k = min(17, n);
searcher = KDTreeSearcher(points);
[neighbors, distances] = knnsearch(searcher, points, 'K', k);
if k < 17
    neighbors(:, end + 1:17) = neighbors(:, end);
    distances(:, end + 1:17) = inf;
end
radii = [1.5, 2.5, 4.0, 6.0, 10.0];
counts = localRadiusCounts(searcher, points, radii, batchSize);
pool.nn1 = min(distances(:, 2), 10.0);
pool.count5 = min(counts(:, 4), 10.0);

center = median(points, 1);
radius = sqrt(sum((points - center).^2, 2));
radiusScale = max(double(prctile(radius, 95)), 1.0);
baseRanks = localAverageRanks(pool.baseScore) ./ max(n, 1);
tissueRanks = localAverageRanks(pool.tissue) ./ max(n, 1);
neighborBase = pool.baseScore(neighbors);
features = single([pool.baseScore, pool.sourceBase, pool.tissue, ...
    baseRanks, tissueRanks, radius ./ radiusScale, ...
    distances(:, [2, 3, 5, 9, 17]), counts, ...
    mean(neighborBase(:, 2:9), 2), min(neighborBase(:, 2:9), [], 2), ...
    mean(neighborBase(:, 2:17), 2), min(neighborBase(:, 2:17), [], 2), ...
    localNeighborhoodShape(points, neighbors(:, 1:8)), ...
    localNeighborhoodShape(points, neighbors(:, 1:16))]);
end

function counts = localRadiusCounts(searcher, points, radii, batchSize)
n = size(points, 1);
counts = zeros(n, numel(radii));
for first = 1:batchSize:n
    last = min(first + batchSize - 1, n);
    [~, distances] = rangesearch(searcher, points(first:last, :), max(radii));
    for localIndex = 1:numel(distances)
        d = distances{localIndex};
        counts(first + localIndex - 1, :) = ...
            max(sum(d(:) <= radii, 1) - 1, 0);
    end
end
end

function shape = localNeighborhoodShape(points, neighborIndices)
n = size(points, 1);
k = size(neighborIndices, 2);
neighbors = reshape(points(neighborIndices(:), :), [n, k, 3]);
centered = neighbors - mean(neighbors, 2);
x = centered(:, :, 1);
y = centered(:, :, 2);
z = centered(:, :, 3);
scale = max(k - 1, 1);
[small, middle, large] = localSymmetricEigenvalues( ...
    sum(x .* x, 2) ./ scale, sum(y .* y, 2) ./ scale, ...
    sum(z .* z, 2) ./ scale, sum(x .* y, 2) ./ scale, ...
    sum(x .* z, 2) ./ scale, sum(y .* z, 2) ./ scale);
normalizer = max(large, 1e-6);
shape = [(large - middle) ./ normalizer, ...
    (middle - small) ./ normalizer, small ./ normalizer];
end

function features = localDenseFeatures(pool, ct, affine, headCenter)
points = pool.points;
n = size(points, 1);
centerValue = double(localSampleWorld(ct, affine, points));
headRadius = sqrt(sum((points - headCenter).^2, 2));
base = [pool.baseScore, pool.sourceBase, pool.nn1, pool.count5, ...
    pool.tissue, headRadius, centerValue];
directions = localShellDirections;
shellMean25 = [];
shellMean35 = [];
radii = [0.8, 1.6, 2.5, 3.5];
for radiusIndex = 1:numel(radii)
    radius = radii(radiusIndex);
    values = zeros(n, size(directions, 1), 'single');
    for directionIndex = 1:size(directions, 1)
        locations = points + radius .* directions(directionIndex, :);
        values(:, directionIndex) = localSampleWorld(ct, affine, locations);
    end
    shellMean = double(mean(values, 2));
    if radius == 2.5
        shellMean25 = shellMean;
    elseif radius == 3.5
        shellMean35 = shellMean;
    end
    shellShape = localShellShape(values, directions);
    base = [base, double(max(values, [], 2)), shellMean, ...
        double(std(values, 1, 2)), double(localRowQuantile(values, 0.25)), ...
        double(localRowQuantile(values, 0.75)), ...
        double(mean(values >= 1.0, 2)), shellShape]; %#ok<AGROW>
end
base = single([base, centerValue - shellMean25, centerValue - shellMean35]);
ranks = zeros(size(base), 'single');
for column = 1:size(base, 2)
    ranks(:, column) = single(localAverageRanks(base(:, column)) ./ max(n, 1));
end
features = [base, ranks];
end

function directions = localShellDirections
directions = zeros(26, 3);
index = 0;
for x = [-1, 0, 1]
    for y = [-1, 0, 1]
        for z = [-1, 0, 1]
            if x == 0 && y == 0 && z == 0
                continue;
            end
            index = index + 1;
            directions(index, :) = [x, y, z] ./ norm([x, y, z]);
        end
    end
end
end

function shape = localShellShape(values, directions)
% NumPy promotes quantiles and the direction-weighted moments to double.
% Preserve that behavior because ranks of very small eigenvalues are useful
% to the frozen classifier.
values = double(values);
floorValue = localRowQuantile(values, 0.35);
weights = max(values - floorValue, 0);
totals = sum(weights, 2);
valid = totals > 1e-8;
weights(valid, :) = weights(valid, :) ./ totals(valid);
weights(~valid, :) = 0;
x = directions(:, 1)';
y = directions(:, 2)';
z = directions(:, 3)';
[small, middle, large] = localSymmetricEigenvalues( ...
    sum(weights .* (x .* x), 2), sum(weights .* (y .* y), 2), ...
    sum(weights .* (z .* z), 2), sum(weights .* (x .* y), 2), ...
    sum(weights .* (x .* z), 2), sum(weights .* (y .* z), 2));
normalizer = small + middle + large;
shape = zeros(size(values, 1), 3);
shape(valid, :) = [small(valid), middle(valid), large(valid)] ./ normalizer(valid);
end

function [small, middle, large] = localSymmetricEigenvalues(a, d, f, b, c, e)
q = (a + d + f) ./ 3;
p1 = b .* b + c .* c + e .* e;
p2 = (a - q).^2 + (d - q).^2 + (f - q).^2 + 2 .* p1;
p = sqrt(max(p2 ./ 6, 0));
regular = p > 1e-14;
r = zeros(size(q));
if any(regular)
    ar = (a(regular) - q(regular)) ./ p(regular);
    dr = (d(regular) - q(regular)) ./ p(regular);
    fr = (f(regular) - q(regular)) ./ p(regular);
    br = b(regular) ./ p(regular);
    cr = c(regular) ./ p(regular);
    er = e(regular) ./ p(regular);
    determinant = ar .* dr .* fr + 2 .* br .* cr .* er - ...
        ar .* er .* er - dr .* cr .* cr - fr .* br .* br;
    r(regular) = max(min(determinant ./ 2, 1), -1);
end
phi = acos(r) ./ 3;
large = q + 2 .* p .* cos(phi);
small = q + 2 .* p .* cos(phi + 2 .* pi ./ 3);
middle = 3 .* q - small - large;
if any(~regular)
    diagonal = sort([a(~regular), d(~regular), f(~regular)], 2);
    small(~regular) = diagonal(:, 1);
    middle(~regular) = diagonal(:, 2);
    large(~regular) = diagonal(:, 3);
end
end

function q = localRowQuantile(values, probability)
sorted = sort(values, 2);
position = (size(values, 2) - 1) .* probability;
lower = floor(position) + 1;
upper = ceil(position) + 1;
fraction = position - floor(position);
q = sorted(:, lower) .* (1 - fraction) + sorted(:, upper) .* fraction;
end

function ranks = localAverageRanks(values)
values = double(values(:));
[~, ~, group] = unique(values, 'sorted');
counts = accumarray(group, 1);
last = cumsum(counts);
first = last - counts + 1;
average = (first + last) ./ 2;
ranks = average(group);
if any(~isfinite(values))
    error('LeG_autoElecs:NonfiniteFeature', ...
        'A feature rank column contains nonfinite values.');
end
end

function values = localSampleWorld(volume, affine, points)
n = size(points, 1);
if n == 0
    values = zeros(0, 1, 'single');
    return;
end
inverse = inv(affine);
voxels = localApplyAffine(inverse, points) + 1;
volumeSize = size(volume);
volumeSize(end + 1:3) = 1;
tolerance = 1e-6;
valid = all(voxels >= 1 - tolerance, 2) & ...
    voxels(:, 1) <= volumeSize(1) + tolerance & ...
    voxels(:, 2) <= volumeSize(2) + tolerance & ...
    voxels(:, 3) <= volumeSize(3) + tolerance;
values = zeros(n, 1, 'single');
if ~any(valid)
    return;
end
v = voxels(valid, :);
v = min(max(v, 1), volumeSize);
lower = floor(v);
upper = min(lower + 1, volumeSize);
fraction = single(v - lower);
x0 = lower(:, 1); x1 = upper(:, 1);
y0 = lower(:, 2); y1 = upper(:, 2);
z0 = lower(:, 3); z1 = upper(:, 3);
fx = fraction(:, 1); fy = fraction(:, 2); fz = fraction(:, 3);
v000 = volume(sub2ind(volumeSize, x0, y0, z0));
v100 = volume(sub2ind(volumeSize, x1, y0, z0));
v010 = volume(sub2ind(volumeSize, x0, y1, z0));
v110 = volume(sub2ind(volumeSize, x1, y1, z0));
v001 = volume(sub2ind(volumeSize, x0, y0, z1));
v101 = volume(sub2ind(volumeSize, x1, y0, z1));
v011 = volume(sub2ind(volumeSize, x0, y1, z1));
v111 = volume(sub2ind(volumeSize, x1, y1, z1));
interpolated = ...
    v000 .* (1 - fx) .* (1 - fy) .* (1 - fz) + ...
    v100 .* fx .* (1 - fy) .* (1 - fz) + ...
    v010 .* (1 - fx) .* fy .* (1 - fz) + ...
    v110 .* fx .* fy .* (1 - fz) + ...
    v001 .* (1 - fx) .* (1 - fy) .* fz + ...
    v101 .* fx .* (1 - fy) .* fz + ...
    v011 .* (1 - fx) .* fy .* fz + ...
    v111 .* fx .* fy .* fz;
values(valid) = single(interpolated);
end

function transformed = localApplyAffine(affine, points)
transformed = [double(points), ones(size(points, 1), 1)] * double(affine)';
transformed = transformed(:, 1:3);
end

function scores = localPredictLightGBM(features, model)
n = size(features, 1);
raw = zeros(n, 1);
isLeaf = logical(model.isLeaf(:));
featureIndex = double(model.featureIndex(:));
threshold = double(model.threshold(:));
leftChild = int32(model.leftChild(:));
rightChild = int32(model.rightChild(:));
defaultLeft = logical(model.defaultLeft(:));
leafValue = double(model.leafValue(:));
roots = int32(model.roots(:));
for tree = 1:numel(roots)
    node = repmat(roots(tree), n, 1);
    while true
        active = find(~isLeaf(double(node)));
        if isempty(active)
            break;
        end
        activeNode = double(node(active));
        columns = featureIndex(activeNode);
        linear = sub2ind(size(features), active, columns);
        values = features(linear);
        goLeft = double(values) <= threshold(activeNode);
        missing = ~isfinite(values);
        goLeft(missing) = defaultLeft(activeNode(missing));
        leftRows = active(goLeft);
        rightRows = active(~goLeft);
        node(leftRows) = leftChild(activeNode(goLeft));
        node(rightRows) = rightChild(activeNode(~goLeft));
    end
    raw = raw + leafValue(double(node));
end
sigmoid = 1;
if isfield(model, 'sigmoid')
    sigmoid = double(model.sigmoid);
end
scores = 1 ./ (1 + exp(-sigmoid .* raw));
end

function selected = localSpatialNms(points, scores, separationMm)
if isempty(points)
    selected = zeros(0, 1);
    return;
end
[~, order] = sort(scores, 'descend');
selected = zeros(numel(order), 1);
count = 0;
for index = 1:numel(order)
    candidate = order(index);
    if count > 0
        distance = sqrt(sum((points(selected(1:count), :) - points(candidate, :)).^2, 2));
        if any(distance < separationMm)
            continue;
        end
    end
    count = count + 1;
    selected(count) = candidate;
end
selected = selected(1:count);
end

function names = localFeatureNames
names = {'base_score'; 'source_base'; 'tissue_probability'; ...
    'rank_base_score'; 'rank_tissue_probability'; 'robust_radius'};
for k = [1, 2, 4, 8, 16]
    names{end + 1, 1} = sprintf('knn_distance_%d', k); %#ok<AGROW>
end
for radius = [1.5, 2.5, 4.0, 6.0, 10.0]
    names{end + 1, 1} = sprintf('neighbor_count_r%gmm', radius); %#ok<AGROW>
end
names = [names; {'neighbor_base_mean_k8'; 'neighbor_base_min_k8'; ...
    'neighbor_base_mean_k16'; 'neighbor_base_min_k16'; ...
    'local_linearity_k8'; 'local_planarity_k8'; 'local_scattering_k8'; ...
    'local_linearity_k16'; 'local_planarity_k16'; 'local_scattering_k16'}];
dense = {'base_score'; 'source_base'; 'nn1_mm'; 'count5'; ...
    'tissue_probability'; 'head_radius_mm'; 'ct_center'};
statistics = {'max', 'mean', 'std', 'q25', 'q75', 'bright_fraction', ...
    'shape_small', 'shape_middle', 'shape_large'};
for radius = [0.8, 1.6, 2.5, 3.5]
    for statistic = 1:numel(statistics)
        dense{end + 1, 1} = sprintf('ct_%s_r%gmm', ...
            statistics{statistic}, radius); %#ok<AGROW>
    end
end
dense = [dense; {'core_shell_contrast_r2p5mm'; 'core_shell_contrast_r3p5mm'}];
dense = [dense; strcat('rank_', dense)];
names = [names; strcat('ct_', dense)];
end

function [WC, label] = localOutputVoxels(app, predictionMm, ctAffine)
mrInfo = localReadApp(app, 'MRInfo', []);
if isstruct(mrInfo) && isfield(mrInfo, 'mat') && isequal(size(mrInfo.mat), [4, 4])
    homogeneous = double(mrInfo.mat) \ [predictionMm, ones(size(predictionMm, 1), 1)]';
    WC = round(homogeneous(1:3, :)');
    label = 'app.MRInfo.mat';
else
    voxelZero = localApplyAffine(inv(ctAffine), predictionMm);
    WC = round(voxelZero + 1);
    label = 'CT NIfTI voxel indices';
end
end

function debug = localEmptyDebug(subject, ctPath, maskPaths, options, model, started)
debug = struct('mode', 'population_contact_classifier_matlab', ...
    'subject', subject, 'final_mm', zeros(0, 3), 'wc', zeros(0, 3), ...
    'selected_scores', zeros(0, 1), 'candidate_mm', zeros(0, 3), ...
    'candidate_scores', zeros(0, 1), 'reference_loaded', false, ...
    'contact_count_used', false, 'patient_identity_used', false, ...
    'ct_path', ctPath, 'mask_paths', {maskPaths}, ...
    'model_path', options.modelPath, ...
    'model_training_scope', localText(model.trainingScope), ...
    'parameters', options, 'total_seconds', toc(started));
end

function localSaveResult(path, predictionMm, selectedScores, debug)
folder = fileparts(path);
if ~isempty(folder) && ~isfolder(folder)
    mkdir(folder);
end
referenceLoaded = false;
contactCountUsed = false;
patientIdentityUsed = false;
save(path, 'predictionMm', 'selectedScores', 'referenceLoaded', ...
    'contactCountUsed', 'patientIdentityUsed', 'debug', '-v7.3');
end

function value = localReadApp(app, name, defaultValue)
value = defaultValue;
try
    if isstruct(app) && isfield(app, name)
        value = app.(name);
    elseif isobject(app) && isprop(app, name)
        value = app.(name);
    end
catch
    value = defaultValue;
end
end

function output = localCellstr(value)
if iscell(value)
    output = cellfun(@localText, value, 'UniformOutput', false);
else
    output = cellstr(string(value));
end
end

function output = localText(value)
if isempty(value)
    output = '';
elseif ischar(value)
    output = value;
elseif isstring(value)
    output = char(value(1));
else
    output = char(string(value));
end
end

function localPrint(options, varargin)
if options.verbose
    fprintf(varargin{:});
end
end
