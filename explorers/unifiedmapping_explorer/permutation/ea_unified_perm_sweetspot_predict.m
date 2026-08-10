function predresults = ea_unified_perm_sweetspot_predict(obj, trainPermtestFile, sigMode, tail, batchSize)
% Out-of-sample validation of a voxelwise permutation test (built via
% ea_unified_perm_sweetspot on a *different*, training ea_unifiedmapping
% object) against this (test) object's own, never-permuted patient scores.
% Each permuted (and the empirical) training map is reprojected onto this
% object's own efield grid via a one-time nearest-neighbor index map -- the
% geometric correspondence only depends on the two objects'
% results.sweetspotmapping.space, not on the map values, so it is computed
% once and reused for every permutation instead of reslicing Nperm+1 times.
%
% trainPermtestFile - path to a results.mat saved by ea_unified_perm_sweetspot.
%           Rperm/pperm are read directly from it via matfile row-chunked
%           reads (it's already saved flat and single-precision -- see
%           ea_unified_perm_sweetspot.m -- so there is no separate
%           conversion step here).
% sigMode - uncorrected significance thresholding (using obj.alphalevel,
%           evaluated on this -- the test -- object) before prediction:
%   'None'        - no thresholding, raw R everywhere (default)
%   'Independent' - each map (empirical & every permutation) is thresholded
%                   using its own p-values
%   'Fixed'       - the empirical map's significant-voxel mask is computed
%                   once and applied to every permutation
% tail    - which direction(s) of the null distribution count as "as
%           extreme as" the empirical prediction, for exceedCount/pperm:
%   'right' (default) - Rpredperm >= Rpredemp (matches a hypothesis that
%                        the map should only ever predict IMPROVEMENT)
%   'left'             - Rpredperm <= Rpredemp
%   'both'             - abs(Rpredperm) >= abs(Rpredemp) (a priori the
%                        projected map could fail to predict in either
%                        direction)
% batchSize - permutations processed per batch (matfile row-chunked read).
%             Default 250; lower it further on tighter-memory machines, at
%             the cost of more, smaller disk reads.

% This feature only supports sweetspotmapping today -- see the matching
% check and rationale in ea_unified_perm_sweetspot.m (obj.drawTool reflects
% GUI display state, not necessarily what the caller means to test; this
% object could have fiberfiltering/networkmapping results too).
if ~isfield(obj.results, 'sweetspotmapping')
    ea_error(['predpermtest currently only supports sweetspotmapping (fiberfiltering/networkmapping ', ...
        'are not yet implemented). obj.results.sweetspotmapping was not found -- compute ', ...
        'sweetspot results for this (test) object first.']);
end

if size(obj.responsevar, 2) > 1
    ea_error('Hemiscore responsevar (2 columns) is not yet supported for permutation testing.');
end

if ~exist('sigMode', 'var') || isempty(sigMode)
    sigMode = 'None';
end
if ~ismember(sigMode, {'None', 'Independent', 'Fixed'})
    ea_error('sigMode must be ''None'', ''Independent'', or ''Fixed''.');
end

if ~exist('tail', 'var') || isempty(tail)
    tail = 'right';
end
if ~ismember(tail, {'right', 'left', 'both'})
    ea_error('tail must be ''right'', ''left'', or ''both''.');
end

if ~exist('batchSize', 'var') || isempty(batchSize)
    batchSize = 250;
end

trainMeta = load(trainPermtestFile, 'space', 'Remp', 'pemp', 'Nperm', 'corrType', 'PermIdx', 'ID', 'nsides');
mTrain = matfile(trainPermtestFile);

nsides = numel(obj.results.sweetspotmapping.space);
srcIdx = cell(1, nsides);
testUncovered = zeros(1, nsides);
testTotal = zeros(1, nsides);
trainUncovered = zeros(1, nsides);
trainTotal = zeros(1, nsides);

fprintf('Voxel coverage report (training grid vs. this object''s test grid):\n');
for side = 1:nsides
    srcIdx{side} = ea_unified_perm_nnindexmap(trainMeta.space{side}, obj.results.sweetspotmapping.space{side});

    testTotal(side) = numel(srcIdx{side});
    testUncovered(side) = sum(isnan(srcIdx{side}));

    trainTotal(side) = numel(trainMeta.space{side}.img);
    covered = unique(srcIdx{side}(~isnan(srcIdx{side})));
    trainUncovered(side) = trainTotal(side) - numel(covered);

    fprintf(['  Side %d: %d/%d (%.1f%%) test-grid voxels have no corresponding training voxel.\n', ...
        '           %d/%d (%.1f%%) training-grid voxels are not represented anywhere in the test grid.\n'], ...
        side, testUncovered(side), testTotal(side), 100*testUncovered(side)/testTotal(side), ...
        trainUncovered(side), trainTotal(side), 100*trainUncovered(side)/trainTotal(side));
end

if isempty(obj.customselection)
    patsel = obj.patientselection;
else
    patsel = obj.customselection;
end
Nperm = trainMeta.Nperm;

% Uncorrected significance thresholding (training-grid space, before
% reprojection -- Remp/pemp and Rperm/pperm share the same voxel indexing,
% so this is a plain elementwise mask).
alphalevel = obj.alphalevel;
Remp_train = trainMeta.Remp;
if ~strcmp(sigMode, 'None')
    for side = 1:nsides
        nonsig = isnan(trainMeta.pemp{side}) | trainMeta.pemp{side} > alphalevel;
        Remp_train{side}(nonsig) = nan;
    end
end
fixedMask = cell(1, nsides); % only used when sigMode == 'Fixed'
if strcmp(sigMode, 'Fixed')
    for side = 1:nsides
        fixedMask{side} = ~isnan(Remp_train{side}); % logical, training-grid space
    end
end

% validMask/efieldT are both independent of which map (empirical or which
% permutation) is being predicted -- precomputed once and reused below
% instead of being recomputed on every one of the Nperm+1 calls.
validMask = cell(1, nsides);
efieldT = cell(1, nsides);
for side = 1:nsides
    validMask{side} = ~isnan(srcIdx{side});
    efieldT{side} = obj.results.sweetspotmapping.efield{side}(patsel,:)';
end

fprintf('Predicting from empirical (unpermuted) training map...\n');
empvals = cell(1, nsides);
for side = 1:nsides
    empvals{side} = nan(numel(srcIdx{side}), 1);
    empvals{side}(validMask{side}) = Remp_train{side}(srcIdx{side}(validMask{side}));
end
[Rpredemp, Rpredemp_pval, Ihatemp] = ea_unified_perm_sweetspot_predictfromvals(obj, empvals, patsel, trainMeta.corrType, efieldT);
if isnan(Rpredemp)
    ea_error('Empirical out-of-sample prediction is NaN -- cannot rank against the null distribution. Check basepredictionon/posvisible/negvisible settings and voxel coverage.');
end
fprintf('Empirical out-of-sample prediction: Rpredemp = %.4f (parametric p = %.4g).\n', Rpredemp, Rpredemp_pval);
fprintf('Check this against your independently-computed empirical R now -- if it does not match, stop here (Ctrl+C) before the %d-permutation null distribution runs.\n', Nperm);

fprintf('Predicting from %d permuted training maps, in batches of %d (serial -- no parfor)...\n', Nperm, batchSize);
Rpredperm = nan(Nperm, 1);
corrType = trainMeta.corrType;

% Serial, not parfor: a persistent worker pool driven by dozens of
% sequential parfor calls (one per batch) has been observed to accumulate
% memory across calls until a worker gets OOM-killed and the whole pool
% fails to recover mid-run -- a different, harder-to-fix failure mode than
% the blanket-load() crash the batching itself solves. Running serially
% removes the worker pool from the picture entirely, at the cost of using
% one core instead of a capped pool.
for batchStart = 1:batchSize:Nperm
    batchIdx = batchStart:min(batchStart+batchSize-1, Nperm);
    nb = numel(batchIdx);

    % Only this batch's rows ever touch memory -- mTrain.Rperm_sideN is a
    % plain top-level array (not nested in a cell), so matfile reads
    % exactly these rows from disk instead of materializing the full
    % Nperm x Nvoxels array.
    Rbatch1 = mTrain.Rperm_side1(batchIdx,:);
    pbatch1 = mTrain.pperm_side1(batchIdx,:);
    if nsides >= 2
        Rbatch2 = mTrain.Rperm_side2(batchIdx,:);
        pbatch2 = mTrain.pperm_side2(batchIdx,:);
    else
        Rbatch2 = [];
        pbatch2 = [];
    end

    RpredBatch = nan(nb, 1);
    for bi = 1:nb
        permvals = cell(1, nsides);
        for side = 1:nsides
            permvals{side} = nan(numel(srcIdx{side}), 1);
            valid = validMask{side};
            if side == 1
                row = double(Rbatch1(bi,:));
                prow = double(pbatch1(bi,:));
            else
                row = double(Rbatch2(bi,:));
                prow = double(pbatch2(bi,:));
            end
            switch sigMode
                case 'Independent'
                    row(isnan(prow) | prow > alphalevel) = nan;
                case 'Fixed'
                    row(~fixedMask{side}) = nan;
            end
            permvals{side}(valid) = row(srcIdx{side}(valid));
        end
        RpredBatch(bi) = ea_unified_perm_sweetspot_predictfromvals(obj, permvals, patsel, corrType, efieldT);
    end
    Rpredperm(batchIdx) = RpredBatch;

    fprintf('  Completed permutations %d-%d of %d\n', batchIdx(1), batchIdx(end), Nperm);
end

% A permutation whose map produced no defined prediction (NaN -- e.g. zero
% significant voxels under sigMode='Independent') is itself evidence
% against the null being able to predict, not missing data: it stays in
% the denominator and counts as not exceeding the empirical result
% (comparisons against NaN are always false in MATLAB, which already gives
% this behavior -- made explicit here rather than left implicit, and
% reported instead of silent).
nNaNperm = sum(isnan(Rpredperm));
switch tail
    case 'right'
        exceedCount = sum(Rpredperm >= Rpredemp);
    case 'left'
        exceedCount = sum(Rpredperm <= Rpredemp);
    case 'both'
        exceedCount = sum(abs(Rpredperm) >= abs(Rpredemp));
end
pperm = exceedCount / Nperm;
fprintf('%d/%d permutations (%.1f%%) produced no defined prediction (NaN) -- counted as not exceeding the empirical result.\n', nNaNperm, Nperm, 100*nNaNperm/Nperm);
disp(['Out-of-sample permuted p (tail=''', tail, ''') = ', sprintf('%0.3f', pperm), ' (empirical R ranks ', num2str(exceedCount), ' of ', num2str(Nperm), ').']);

predresults.Rpredemp = Rpredemp;
predresults.Rpredemp_pval = Rpredemp_pval; % parametric p-value of the empirical correlation itself (distinct from pperm, the permutation-based null p-value)
predresults.Ihatemp = Ihatemp; % per-patient predicted score underlying Rpredemp
predresults.Iemp = obj.responsevar(patsel); % this object's real, unpermuted scores -- both saved so the correlation plot can be regenerated standalone later, without a live object
predresults.responsevarlabel = obj.responsevarlabel;
predresults.Rpredperm = Rpredperm;
predresults.nNaNperm = nNaNperm; % permutations with no defined prediction, counted as not exceeding Rpredemp (see pperm)
predresults.exceedCount = exceedCount; % how many (of Nperm) permutations were as extreme as Rpredemp per `tail`; pperm = exceedCount/Nperm
predresults.pperm = pperm;
predresults.tail = tail;
predresults.trainPermIdx = trainMeta.PermIdx; % training cohort's shuffle indices (NOT the test cohort -- this object's own patients are never permuted). Columns correspond 1:1 to Rpredperm entries, traceable to the training permtest that produced each permuted map.
predresults.trainPermtestFile = trainPermtestFile;
predresults.trainID = trainMeta.ID;
predresults.testID = obj.ID;
predresults.corrType = trainMeta.corrType;
predresults.basepredictionon = obj.basepredictionon;
predresults.sigMode = sigMode;
predresults.alphalevel = alphalevel;
predresults.testUncoveredVoxels = testUncovered;
predresults.testTotalVoxels = testTotal;
predresults.trainUncoveredVoxels = trainUncovered;
predresults.trainTotalVoxels = trainTotal;

predresults.settings.sigMode = sigMode;
predresults.settings.tail = tail;
predresults.settings.batchSize = batchSize;
predresults.settings.alphalevel = alphalevel;
predresults.settings.basepredictionon = obj.basepredictionon;
predresults.settings.posvisible = obj.posvisible;
predresults.settings.negvisible = obj.negvisible;

% Filename identifies which training run this test cohort was checked
% against (a test cohort may be validated against several different
% training analyses for different purposes -- each gets its own folder
% instead of overwriting the last one), plus sigMode and a timestamp. Own
% self-contained folder, same convention as ea_unified_perm_sweetspot.m.
[~, trainRunName] = fileparts(fileparts(trainPermtestFile)); % trainPermtestFile lives in its own run folder; use that folder's name
if isempty(trainRunName)
    [~, trainRunName] = fileparts(trainPermtestFile);
end
runName = sprintf('%s_predpermtest_vs_%s_%s_%s', obj.ID, trainRunName, sigMode, datestr(now, 'yyyymmdd_HHMMSS'));
outdir = fullfile(fileparts(obj.leadgroup), 'UnifiedMappingExplorer', 'permutation', runName, filesep);
suffix = 1;
while exist(outdir, 'dir') % timestamp makes this vanishingly rare, but never silently overwrite/merge into an existing run's folder
    suffix = suffix + 1;
    outdir = fullfile(fileparts(obj.leadgroup), 'UnifiedMappingExplorer', 'permutation', sprintf('%s_%d', runName, suffix), filesep);
end
ea_mkdir(outdir);

outfile = fullfile(outdir, 'results.mat');
predresults.savedfile = outfile;
save(outfile, '-struct', 'predresults', '-v7.3');
fprintf('Saved out-of-sample permutation prediction results to %s\n', outfile);

ea_unified_perm_sweetspot_predict_plot(outfile); % shared with standalone re-plotting from a saved file
