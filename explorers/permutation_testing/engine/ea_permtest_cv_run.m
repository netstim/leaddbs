function run = ea_permtest_cv_run(cvb, permset, opts)
% Prediction permutation test ("Q5"): is the explorer's k-fold
% cross-validated prediction better than chance? The whole k-fold pipeline
% (fiber model per training fold incl. display threshold -> prediction of
% the test patients, see ea_permtest_cv_predict.m) is repeated for every
% shuffle of the scores, with the same fold splits as the real run.
%
% Reported, each with Spearman and Pearson R of observed vs predicted
% (complete cases, as the explorer's CV plot):
%   (b) every model against its own score (subscore, or the response
%       variable), one-sided p = (#shuffles with R >= real R + 1) / (Nperm + 1),
%       uncorrected (each subscore is its own hypothesis); the subscores
%       are shuffled jointly, a patient's whole row moves together
%   (a) optional: the combined (weighted) prediction of all subscore models
%       against a chosen clinical variable, which is shuffled together with
%       the subscores (as are patient-specific weight variables)
% The explorer's 'Post-hoc correct for group' is applied to every
% prediction when requested (and the lead group has several cohorts).
%
% cvb     - from ea_permtest_cv_adapter_disctract (checked with
%           ea_permtest_cv_check_disctract).
% permset - from ea_permtest_get_permset, matching cvb.patientIDs.
% opts    - struct, all optional:
%   .name            analysis name (see ea_permtest_setup.m). Default analysis_<time>.
%   .stratifyByGroup shuffle only within cohort groups. Default false.
%   .batchSize       shuffles computed together. Default 50.
%   .confirmFcn      function handle @(message) -> true/false, asked with
%                    the estimated duration before the run starts when
%                    k-fold iterations > 1 (e.g. a dialog). Default: no question.
%   .checkReport     check report (documented in the README).
%   .permsetFile     permset file (documented in the README).
%
% Run folder: <cvb.resultsdir>/<name>/prediction_k<k>[_iter<n>]/
%   README.md, run_info.mat (settings, folds, shuffle info, check, timings),
%   results.mat (predictions, R and p per shuffle), null-distribution figures.
%
% run - struct: .rundir, .cancelled, .results (see results.mat), .seconds

if ~exist('opts', 'var') || isempty(opts)
    opts = struct;
end
opts = local_defaults(opts, struct('name', '', 'stratifyByGroup', false, 'batchSize', 50, ...
    'confirmFcn', [], 'checkReport', [], 'permsetFile', ''));
if isempty(opts.name)
    opts.name = ['analysis_', datestr(now, 'yyyymmdd_HHMMSS')];
    fprintf('No analysis name given, using ''%s''.\n', opts.name);
end
if ~isequal(permset.patientIDs, cvb.patientIDs)
    ea_error('The permset does not belong to this lead group (patient lists differ).');
end
if isstruct(opts.checkReport) && isfield(opts.checkReport, 'ok') && ~opts.checkReport.ok
    ea_error('The check against the explorer''s cross-validation failed -- no permutations were calculated.');
end

cv = cvb.cv;
V = numel(cvb.voterLabels);
sel = cvb.patientselection(:)';
Nsel = numel(sel);
hasTarget = ~isempty(cvb.extra.target);
run = struct('rundir', '', 'cancelled', false, 'results', [], 'seconds', NaN);

% Shuffles: jointly over the scores in use (models, target, weight variables)
used = 1:V;
if hasTarget
    used(end+1) = V + 1;
end
used = [used, V + 1 + find(isnan(cv.weights))];
O0 = ea_permtest_cv_scores(cvb);
strata = [];
if opts.stratifyByGroup
    strata = cvb.group;
end
[~, PermIdxUsed, shuffleInfo] = ea_permtest_derive_shuffles(permset, O0(:, used), sel, strata);
Nperm = size(PermIdxUsed, 2);
PermIdx = repmat((1:cvb.Npatients)', 1, Nperm, size(O0, 2));
PermIdx(:, :, used) = reshape(PermIdxUsed, cvb.Npatients, Nperm, numel(used));

% Timing probe: real scores + a few shuffles, first k-fold iteration
nProbe = min(4, Nperm);
tProbe = tic;
ea_permtest_cv_predict(cvb, cat(3, O0, ea_permtest_cv_scores(cvb, PermIdx(:, 1:nProbe, :))), 1);
perVersion = toc(tProbe) / (nProbe + 1);
estimate = perVersion * (Nperm + 1) * cv.kIter;
msg = sprintf(['Estimated duration: %s (%d shuffles x %d k-fold iteration(s), k = %d, %d model(s)).\n', ...
    'MATLAB is busy until the run has finished.'], local_duration(estimate), Nperm, cv.kIter, cv.kfold, V);
fprintf('%s\n', msg);
if cv.kIter > 1 && ~isempty(opts.confirmFcn)
    if ~opts.confirmFcn(sprintf(['k-fold iterations > 1 repeat the whole cross-validation for every iteration ', ...
            'and every shuffle -- this may take long.\n\n%s\n\nContinue?'], msg))
        run.cancelled = true;
        fprintf('Cancelled -- nothing was saved.\n');
        return
    end
end

% Run folder
runFolder = sprintf('prediction_k%d', cv.kfold);
if cv.kIter > 1
    runFolder = sprintf('%s_iter%d', runFolder, cv.kIter);
end
[rundir, setupdir] = ea_permtest_setup(cvb.resultsdir, opts.name, cvb, permset, opts.stratifyByGroup, runFolder);
fprintf('Prediction permutation run folder:\n%s\n', rundir);

% Real + shuffled predictions, in batches (version 1 = real scores)
tStart = tic;
voterPred = nan(Nsel, V, Nperm + 1);
combPred = nan(Nsel, Nperm + 1);
nModel = zeros(cv.kIter, cv.kfold, V, cvb.Nsides, Nperm + 1, 'uint32');
nAnal = [];
allFast = true;
for b0 = 0:opts.batchSize:Nperm
    idx = b0:min(b0 + opts.batchSize - 1, Nperm); % 0 = real scores
    Ob = nan(size(O0, 1), size(O0, 2), numel(idx));
    if idx(1) == 0
        Ob(:, :, 1) = O0;
    end
    p = idx(idx > 0);
    if ~isempty(p)
        Ob(:, :, idx > 0) = ea_permtest_cv_scores(cvb, PermIdx(:, p, :));
    end
    out = ea_permtest_cv_predict(cvb, Ob);
    voterPred(:, :, idx + 1) = out.voterPred;
    combPred(:, idx + 1) = out.combPred;
    nModel(:, :, :, :, idx + 1) = out.nModelFibers;
    nAnal = out.nAnalysed;
    allFast = allFast && out.fast;
    elapsed = toc(tStart);
    fprintf('  %d/%d shuffles (%.0f s, ~%.0f s left)\n', idx(end), Nperm, elapsed, ...
        elapsed / (idx(end) + 1) * (Nperm - idx(end)));
end
seconds = toc(tStart);

% Post-hoc correction for cohort (explorer option)
if cv.posthoccorrect
    reg = ea_cohortregressor(cvb.group(sel));
    for k = 1:Nperm + 1
        combPred(:, k) = ea_resid(reg, combPred(:, k));
        for v = 1:V
            voterPred(:, v, k) = ea_resid(reg, voterPred(:, v, k));
        end
    end
end

% R of observed vs predicted, for every version (observed = that version's scores)
types = {'Spearman', 'Pearson'};
Rvoter = nan(V, 2, Nperm + 1);
Rcomb = nan(2, Nperm + 1);
for k = 1:Nperm + 1
    if k == 1
        Ok = O0;
    else
        Ok = ea_permtest_cv_scores(cvb, PermIdx(:, k - 1, :));
    end
    for t = 1:2
        for v = 1:V
            Rvoter(v, t, k) = corr(Ok(sel, v), voterPred(:, v, k), 'rows', 'complete', 'type', types{t});
        end
        if hasTarget
            Rcomb(t, k) = corr(Ok(sel, V + 1), combPred(:, k), 'rows', 'complete', 'type', types{t});
        end
    end
end

% p-values: one-sided (prediction better than chance)
res = struct;
res.patients = cvb.patientIDs(sel);
res.voterLabels = cvb.voterLabels;
res.targetLabel = cvb.extra.targetLabel;
res.types = types;
res.folds = cv.folds;
res.observed = O0(sel, 1:V);
res.voterPredReal = voterPred(:, :, 1);
res.voterPredPerm = single(voterPred(:, :, 2:end));
res.combPredReal = combPred(:, 1);
res.combPredPerm = single(combPred(:, 2:end));
res.Rvoter = Rvoter(:, :, 1);
res.RvoterPerm = reshape(Rvoter(:, :, 2:end), V, 2, Nperm);
res.pvoter = nan(V, 2);
for t = 1:2
    for v = 1:V
        res.pvoter(v, t) = local_p(squeeze(Rvoter(v, t, 2:end)), Rvoter(v, t, 1));
    end
end
if hasTarget
    res.observedTarget = O0(sel, V + 1);
    res.Rcomb = Rcomb(:, 1);
    res.RcombPerm = Rcomb(:, 2:end);
    res.pcomb = [local_p(Rcomb(1, 2:end), Rcomb(1, 1)); local_p(Rcomb(2, 2:end), Rcomb(2, 1))];
end
res.nModelFibers = nModel;
res.nAnalysed = nAnal;
save(fullfile(rundir, 'results.mat'), '-struct', 'res', '-v7.3');

settings = rmfield(cvb, 'fib'); %#ok<NASGU>
runOptions = rmfield(opts, {'checkReport', 'confirmFcn'}); %#ok<NASGU>
checkReport = opts.checkReport; %#ok<NASGU>
permsetFingerprint = permset.fingerprint; %#ok<NASGU>
permsetFile = opts.permsetFile; %#ok<NASGU>
fast = allFast; %#ok<NASGU>
save(fullfile(rundir, 'run_info.mat'), 'settings', 'runOptions', 'checkReport', 'shuffleInfo', ...
    'permsetFingerprint', 'permsetFile', 'seconds', 'fast', '-v7.3');

local_readme(rundir, cvb, permset, opts, shuffleInfo, res, seconds, allFast, hasTarget);
local_figures(rundir, res, hasTarget);
fprintf('\nDone in %.0f s. Results in:\n%s\n', seconds, rundir);

run.rundir = rundir;
run.setupdir = setupdir;
run.results = res;
run.seconds = seconds;


function p = local_p(null, real)
% (#shuffles >= real + 1) / (Nperm + 1); NaN never exceeds but counts
if isnan(real)
    p = NaN;
    return
end
p = (sum(null(:) >= real) + 1) / (numel(null) + 1);


function s = local_duration(sec)
if sec < 120
    s = sprintf('%.0f s', sec);
elseif sec < 7200
    s = sprintf('%.0f min', sec / 60);
else
    s = sprintf('%.1f h', sec / 3600);
end


function local_readme(rundir, cvb, permset, opts, info, res, seconds, allFast, hasTarget)
cv = cvb.cv;
L = {};
L{end+1} = sprintf('- Analysis: `%s`, explorer ID: `%s`, mode: %s', opts.name, cvb.ID, cvb.mode);
L{end+1} = sprintf('- Cross-validation: k-fold (randomized), k = %d, iterations: %d, fold splits from the explorer seed (`%s`)', ...
    cv.kfold, cv.kIter, char(string(cvb.rngseed)));
L{end+1} = sprintf('- Patients: %d selected, mirrored E-fields: %d', numel(cvb.patientselection), cvb.mirrorsides);
L{end+1} = sprintf('- Fiber test: %s, coverage threshold %g%%, E-field threshold %g', cvb.stattest, ...
    cvb.statsettings.connthreshold, cvb.statsettings.efieldthreshold);
if cv.showsignificantonly
    L{end+1} = sprintf('- Significance masking: on (%s, alpha %g)', cv.multcompstrategy, cv.alphalevel);
else
    L{end+1} = '- Significance masking: off';
end
L{end+1} = sprintf('- Display threshold: %s; per model (pos/neg visible, amount right/left): %s', cv.threshstrategy, ...
    strjoin(arrayfun(@(v) sprintf('%s %d/%d %g/%g', cvb.voterLabels{v}, cv.vis.pos(v), cv.vis.neg(v), ...
    cv.vis.posamount(v, 1), cv.vis.posamount(v, 2)), 1:numel(cvb.voterLabels), 'Uni', 0), '; '));
L{end+1} = sprintf('- Prediction: %s, model normalization: %s, fit to scores: %d (%s)', cv.basepredictionon, ...
    cv.modelNormalization, cv.doactualprediction, cv.predictionmodel);
if cv.posthoccorrect
    L{end+1} = '- Post-hoc correction for cohort: on';
elseif cv.posthoccorrectRequested
    L{end+1} = '- Post-hoc correction for cohort: requested, not applied (only one cohort in the lead group, as in the explorer)';
else
    L{end+1} = '- Post-hoc correction for cohort: off';
end
if cv.subscoreMode
    L{end+1} = sprintf('- Subscore weights: %s', strjoin(arrayfun(@(w) local_weight(w), cv.weights, 'Uni', 0), ', '));
end
if opts.stratifyByGroup
    L{end+1} = sprintf('- Shuffles: within cohort groups (stratified), %d permutations, permset `%s`, jointness %.2f', ...
        permset.Nperm, permset.fingerprint, info.jointness);
else
    L{end+1} = sprintf('- Shuffles: across all selected patients (pooled), %d permutations, permset `%s`, jointness %.2f', ...
        permset.Nperm, permset.fingerprint, info.jointness);
end
if isstruct(opts.checkReport) && isfield(opts.checkReport, 'ok')
    if opts.checkReport.ok
        L{end+1} = sprintf('- Check against the explorer''s own k-fold CV: **passed** (%d patients, max |diff| %.1e)', ...
            opts.checkReport.Ncompared, opts.checkReport.maxdiff);
    else
        L{end+1} = '- Check against the explorer''s own k-fold CV: **FAILED**';
    end
end
if allFast
    comp = 'fast path for the fiber tests';
else
    comp = 'test function for (some) fiber tests';
end
L{end+1} = sprintf('- Computation: %.0f s, %s', seconds, comp);
L{end+1} = '';
L{end+1} = '| Prediction | Observed | Spearman R | p | Pearson R | p |';
L{end+1} = '|---|---|---|---|---|---|';
for v = 1:numel(res.voterLabels)
    L{end+1} = sprintf('| model %s | %s | %.3f | %.4g | %.3f | %.4g |', res.voterLabels{v}, ...
        res.voterLabels{v}, res.Rvoter(v, 1), res.pvoter(v, 1), res.Rvoter(v, 2), res.pvoter(v, 2)); %#ok<AGROW>
end
if hasTarget
    L{end+1} = sprintf('| combined (weighted) | %s | %.3f | %.4g | %.3f | %.4g |', res.targetLabel, ...
        res.Rcomb(1), res.pcomb(1), res.Rcomb(2), res.pcomb(2));
end
L{end+1} = '';
L{end+1} = ['- p: one-sided (prediction better than chance), (#shuffles with R >= real R + 1) / (Nperm + 1), ', ...
    'uncorrected (not corrected for the number of subscores).'];
ea_permtest_readme_append(rundir, 'Prediction permutation test (k-fold cross-validation)', L);


function s = local_weight(w)
if isnan(w)
    s = 'from variable';
else
    s = sprintf('%g', w);
end


function local_figures(rundir, res, hasTarget)
V = numel(res.voterLabels);
for t = 1:2
    n = V + hasTarget;
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 420 * min(n, 3), 320 * ceil(n / 3)]);
    tl = tiledlayout(fig, ceil(n / 3), min(n, 3), 'TileSpacing', 'compact', 'Padding', 'compact');
    for v = 1:V
        ax = nexttile(tl);
        ea_permtest_nulldist_panel(ax, squeeze(res.RvoterPerm(v, t, :)), res.Rvoter(v, t), ...
            sprintf('%s R (observed vs predicted)', res.types{t}), ['Model ', res.voterLabels{v}], false);
    end
    if hasTarget
        ax = nexttile(tl);
        ea_permtest_nulldist_panel(ax, res.RcombPerm(t, :), res.Rcomb(t), ...
            sprintf('%s R (observed vs predicted)', res.types{t}), ['Combined vs ', res.targetLabel], false);
    end
    title(tl, sprintf('k-fold prediction vs shuffled scores (%s)', res.types{t}));
    exportgraphics(fig, fullfile(rundir, sprintf('prediction_null_%s.png', lower(res.types{t}))), 'Resolution', 300);
    close(fig);
end


function opts = local_defaults(opts, defaults)
fn = fieldnames(defaults);
for i = 1:numel(fn)
    if ~isfield(opts, fn{i}) || (isempty(opts.(fn{i})) && ~isstruct(defaults.(fn{i})))
        opts.(fn{i}) = defaults.(fn{i});
    end
end
