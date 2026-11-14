function cvb = ea_permtest_cv_adapter_disctract(obj, pathopts, cvopts)
% Adapter for the prediction permutation test (k-fold cross-validation) of
% the fiber filtering explorer (ea_disctract). Collects everything
% ea_permtest_cv_predict.m needs to repeat the explorer's k-fold CV
% (ea_disctract.kfoldcv -> crossval -> ea_compute_fibscore_model ->
% ea_discfibers_calcstats) for shuffled scores, without touching obj.
%
% Every fold builds the explorer's fiber model from its training patients
% (coverage filter, fiber test, optional significance masking, display
% threshold, optional model normalization) and scores the test patients
% with the explorer's prediction method ("profile of scores", mean, sum,
% peak). Fold splits are the explorer's own (rng(obj.rngseed), cvpartition).
%
% obj     - ea_disctract object (not modified).
% pathopts - see ea_permtest_adapter_disctract.m (.fibfiltPath / .groupdir).
% cvopts  - struct, all optional:
%   .kIter          k-fold iterations (the explorer's 'Iterations' field).
%                   Default obj.kIter.
%   .posthoccorrect correct the predictions for cohort, as the explorer's
%                   'Post-hoc correct for group' option. Default false.
%   .targetLabel    clinical variable (obj.M.clinical.labels) for the
%                   optional combined prediction ('Split & Color By
%                   Subscore' only). Default '' (none).
%
% cvb - the fiber bundle of ea_permtest_adapter_disctract.m (without its
%       fiber sets) plus:
%   .voterLabels   one per model (subscore, or the response variable)
%   .extra         .target (Npatients x 1, [] if none), .targetLabel,
%                  .weightvars (Npatients x Nvoters, NaN columns for fixed weights)
%   .cv            k-fold settings (see below)
%   .fib           .model{side} fiber values for the fiber test (all patient
%                  columns, incl. mirrored), .pred{side} values used to score
%                  patients (columns = patients), .rowidx{side} rows of the
%                  explorer's fibsval (only fibers touched by a selected patient)

if ~exist('pathopts', 'var') || isempty(pathopts)
    pathopts = struct;
end
if ~exist('cvopts', 'var') || isempty(cvopts)
    cvopts = struct;
end
if ~isfield(cvopts, 'kIter') || isempty(cvopts.kIter)
    cvopts.kIter = obj.kIter;
end
if ~isfield(cvopts, 'posthoccorrect') || isempty(cvopts.posthoccorrect)
    cvopts.posthoccorrect = false;
end
if ~isfield(cvopts, 'targetLabel')
    cvopts.targetLabel = '';
end

if ~ismember(obj.multitractmode, {'Single Tract Analysis', 'Split & Color By Subscore'})
    ea_error(sprintf(['The prediction permutation test supports ''Single Tract Analysis'' and ', ...
        '''Split & Color By Subscore'' (current mode: ''%s'').'], obj.multitractmode), 'title', 'Permutation testing');
end
if ~ismember(obj.statsettings.stimulationmodel, {'Electric Field', 'Sigmoid Field'})
    ea_error(['The prediction permutation test supports E-field / sigmoid-field models only (current: ', ...
        obj.statsettings.stimulationmodel, ').'], 'title', 'Permutation testing');
end
if isstruct(obj.ADJ)
    ea_error('Fiber adjacency smoothing (ADJ) is not supported by the prediction permutation test.', 'title', 'Permutation testing');
end
methods = {'profile of scores: spearman', 'profile of scores: pearson', 'profile of scores: bend', ...
    'mean of scores', 'sum of scores', 'peak of scores', 'peak 5% of scores'};
if ~ismember(lower(obj.basepredictionon), methods)
    ea_error(['Prediction method not supported: ', obj.basepredictionon], 'title', 'Permutation testing');
end
if ~ismember(obj.modelNormalization, {'None', 'z-score', 'van Albada 2007'})
    ea_error(['Model normalization not supported: ', obj.modelNormalization], 'title', 'Permutation testing');
end
thresholds = {'Fixed Amount', 'Percentage Relative to Peak', 'Percentage Relative to Amount', ...
    'Histogram (CDF)', 'Fixed Fiber Value'};
if ~ismember(obj.threshstrategy, thresholds)
    ea_error(['Threshold strategy not supported: ', obj.threshstrategy], 'title', 'Permutation testing');
end
if ~isnumeric(obj.kfold) || obj.kfold < 2
    ea_error('k-fold cross-validation needs k >= 2 (k = 1 is the circular, in-sample prediction).', 'title', 'Permutation testing');
end
if cvopts.kIter < 1
    ea_error('The number of k-fold iterations must be at least 1.', 'title', 'Permutation testing');
end

% Fiber bundle (refuses covariates, hemiscores, PCA, unsupported tests)
cvb = ea_permtest_adapter_disctract(obj, pathopts);
cvb = rmfield(cvb, 'fibersets');
subscoreMode = strcmp(obj.multitractmode, 'Split & Color By Subscore');
V = size(cvb.outcomes, 2);
cvb.voterLabels = cvb.outcomelabels;
N = cvb.Npatients;
sel = cvb.patientselection;
Nsel = numel(sel);
if Nsel < obj.kfold
    ea_error(sprintf('Only %d selected patients for %d folds.', Nsel, obj.kfold), 'title', 'Permutation testing');
end

% Combined prediction ('Split & Color By Subscore'): target and weights
extra = struct('target', [], 'targetLabel', '', 'weightvars', nan(N, V));
weights = ones(1, V);
if subscoreMode
    if isempty(obj.subscore.weights) || numel(obj.subscore.weights) < V
        ea_error('The subscore weights are not set in the explorer.', 'title', 'Permutation testing');
    end
    weights = obj.subscore.weights(1:V);
    weights = weights(:)';
    for v = find(isnan(weights))
        w = obj.subscore.weightvars{v};
        if numel(w) ~= N || size(w, 2) > 1
            ea_error(sprintf('Weight variable of subscore %d does not have one value per patient.', v), 'title', 'Permutation testing');
        end
        extra.weightvars(:, v) = w(:);
    end
end
if ~isempty(cvopts.targetLabel)
    if ~subscoreMode
        ea_error('A combined prediction target is only available in ''Split & Color By Subscore'' mode.', 'title', 'Permutation testing');
    end
    [is, ix] = ismember(cvopts.targetLabel, obj.M.clinical.labels);
    if ~is
        ea_error(['Clinical variable not found: ', cvopts.targetLabel], 'title', 'Permutation testing');
    end
    t = obj.M.clinical.vars{ix};
    if numel(t) ~= N || size(t, 2) > 1
        ea_error('The combined target must have one value per patient (no hemiscores).', 'title', 'Permutation testing');
    end
    extra.target = t(:);
    extra.targetLabel = char(cvopts.targetLabel);
end
cvb.extra = extra;

% Fold splits, exactly as ea_disctract.kfoldcv draws them
stream = rng;
rng(obj.rngseed);
folds = zeros(Nsel, cvopts.kIter);
for it = 1:cvopts.kIter
    cvp = cvpartition(Nsel, 'KFold', obj.kfold);
    for c = 1:cvp.NumTestSets
        folds(cvp.test(c), it) = c;
    end
end
rng(stream);

% Display thresholds per model and side, as ea_discfibers_calcstats reads them
vis = struct('pos', false(1, V), 'neg', false(1, V), 'posamount', zeros(V, 2), 'negamount', zeros(V, 2));
for v = 1:V
    if subscoreMode && ~obj.subscore.special_case
        vis.pos(v) = logical(obj.subscore.posvisible(v));
        vis.neg(v) = logical(obj.subscore.negvisible(v));
        vis.posamount(v, :) = obj.subscore.vis.showposamount(v, 1:2);
        vis.negamount(v, :) = obj.subscore.vis.shownegamount(v, 1:2);
    else
        vis.pos(v) = logical(obj.posvisible);
        vis.neg(v) = logical(obj.negvisible);
        vis.posamount(v, :) = obj.showposamount(1:2);
        vis.negamount(v, :) = obj.shownegamount(1:2);
    end
end

groupsAll = cvb.group;
cv = struct;
cv.strategy = 'k-fold (randomized)';
cv.kfold = obj.kfold;
cv.kIter = cvopts.kIter;
cv.folds = folds;
cv.basepredictionon = lower(obj.basepredictionon);
cv.modelNormalization = obj.modelNormalization;
cv.showsignificantonly = logical(obj.showsignificantonly);
cv.multcompstrategy = obj.multcompstrategy;
cv.alphalevel = obj.alphalevel;
cv.threshstrategy = obj.threshstrategy;
cv.vis = vis;
cv.peakvis = [logical(obj.posvisible), logical(obj.negvisible)]; % ea_discfibers_getpeak uses the global switches
cv.doactualprediction = logical(obj.doactualprediction);
cv.predictionmodel = lower(obj.predictionmodel);
cv.weights = weights;
cv.subscoreMode = subscoreMode;
% the explorer only corrects when the lead group has more than one cohort
cv.posthoccorrectRequested = logical(cvopts.posthoccorrect);
cv.posthoccorrect = logical(cvopts.posthoccorrect) && numel(unique(groupsAll)) > 1;
cvb.cv = cv;

% Fiber values: model (as ea_discfibers_calcstats) and prediction (as
% ea_disctract.crossval), restricted to fibers touched by a selected patient
connid = cvb.connid;
methodid = cvb.methodid;
switch obj.statsettings.stimulationmodel
    case 'Sigmoid Field'
        if obj.connectivity_type == 2
            model = cellfun(@full, obj.results.(connid).('PAM_probA').fibsval, 'Uni', 0);
        else
            model = cellfun(@full, obj.results.(connid).(methodid).fibsval, 'Uni', 0);
            for side = 1:numel(model)
                model{side} = ea_SigmoidFromEfield(model{side});
            end
        end
    otherwise
        model = cellfun(@full, obj.results.(connid).(methodid).fibsval, 'Uni', 0);
end
pred = cellfun(@full, obj.results.(connid).(methodid).fibsval, 'Uni', 0);
gcols = sel;
if cvb.mirrorsides
    gcols = [sel, sel + N];
end
fib = struct('model', {cell(1, cvb.Nsides)}, 'pred', {cell(1, cvb.Nsides)}, 'rowidx', {cell(1, cvb.Nsides)});
for side = 1:cvb.Nsides
    rows = find(any(model{side}(:, gcols) ~= 0 & ~isnan(model{side}(:, gcols)), 2));
    fib.rowidx{side} = rows;
    fib.model{side} = model{side}(rows, :);
    fib.pred{side} = pred{side}(rows, 1:N);
end
cvb.fib = fib;
