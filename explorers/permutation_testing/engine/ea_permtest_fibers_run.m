function run = ea_permtest_fibers_run(bundle, permset, opts)
% Fiber permutation test: for every analysis of an adapter bundle (single
% tract, each subscore, or each group) and every side, compute the real
% fiber map and one fiber map per shuffle of the outcome, with the
% explorer's own statistical test (bundle.statfile), and save everything
% to a new, self-documenting run folder. Tool-independent: works on the
% bundle only (see ea_permtest_adapter_*.m).
%
% Shuffles come from the cohort-level permset (ea_permtest_bundle_shuffles):
% subscores of one selection are read from the same chains, groups get
% their own, and the analysed fibers are identical in every permutation.
%
% Spearman and Pearson use a fast path (ea_permtest_fastcorr*.m): the fiber
% values are ranked once and all shuffles of a batch are computed with one
% matrix product, reproducing corr's R and p. Before it is used, it is
% checked against the test function on the real outcome and the first two
% shuffles of every analysis and side; other tests, or a failed check, use
% the test function (in parallel).
%
% bundle  - from an adapter (consistency-checked, see ea_permtest_check_*).
% permset - from ea_permtest_get_permset, matching bundle.patientIDs.
% opts    - struct, all optional:
%   .name             name of the analysis (folder below bundle.resultsdir,
%                     see ea_permtest_setup.m). Default: analysis_<time>.
%   .stratifyByGroup  shuffle only within bundle.group strata. Default false.
%   .fastPath         use the fast path where possible. Default true.
%   .maxWorkers       'auto' (default) or a number, see ea_permtest_parpool
%                     (only used without the fast path).
%   .batchSize        permutations computed before writing to disk. Default 100.
%   .checkReport      consistency-check report, documented in the README.
%   .permsetFile      path of the permset, documented in the README.
%
% Run folder: <bundle.resultsdir>/<name>/fibers/
%   README.md                 settings, shuffles, check, timings
%   run_info.mat              settings (bundle without fiber data), shuffle
%                             info per analysis, permset fingerprint/file,
%                             check report, computation path, timings
%   fibers.mat                analysed fibers per fiber set: explorer row
%                             indices (nonemptyidx), Nfibers, coordinates
%                             (fibcell)
%   results_<label>.mat (-v7.3), one per analysis:
%       label, outcomeCol, outcome (Npatients x 1), patientselection,
%       gpatsel, fiberset, Nperm, PermIdx (Npatients x Nperm),
%       shuffleInfo, statfile, and per side s:
%       fibIdx_side<s>  analysed fibers (rows of the explorer's fibsval{s})
%       Remp_side<s>    real statistic per fiber (Nfib x 1, double)
%       pemp_side<s>    its p-value (Nfib x 1, double)
%       Rperm_side<s>   statistic per permutation and fiber (Nperm x Nfib, single)
%       pperm_side<s>   p-value per permutation and fiber (Nperm x Nfib, single)
%   Rperm/pperm are written in batches and meant to be read in chunks
%   (matfile), so memory use stays small regardless of the size.
%
% run - struct: .rundir, .setupdir, .name, .resultFiles, .nWorkers,
%       .seconds (per analysis), .compute (path used per analysis).

if ~exist('opts', 'var') || isempty(opts)
    opts = struct;
end
opts = local_defaults(opts, struct('name', '', 'stratifyByGroup', false, 'fastPath', true, ...
    'maxWorkers', 'auto', 'batchSize', 100, 'checkReport', [], 'permsetFile', ''));
if isempty(opts.name)
    opts.name = ['analysis_', datestr(now, 'yyyymmdd_HHMMSS')];
    fprintf('No analysis name given, using ''%s''.\n', opts.name);
end

if isfield(bundle, 'covars') && ~isempty(bundle.covars)
    ea_error('Permutation testing is not designed for covariates yet. No permutations were calculated.', ...
        'title', 'Permutation testing');
end
if ~isequal(permset.patientIDs, bundle.patientIDs)
    ea_error('The permset does not belong to this lead group (patient lists differ).');
end
if isstruct(opts.checkReport) && isfield(opts.checkReport, 'ok') && ~opts.checkReport.ok
    ea_error('The consistency check failed -- not running permutations for this analysis.');
end

tStart = tic;
shuffles = ea_permtest_bundle_shuffles(permset, bundle, opts.stratifyByGroup);

% Run folder <resultsdir>/<name>/fibers (never overwrites an earlier run)
[rundir, setupdir] = ea_permtest_setup(bundle.resultsdir, opts.name, bundle, permset, opts.stratifyByGroup, 'fibers');
fprintf('Permutation run folder:\n%s\n', rundir);

% Workers are only started when a test has to run without the fast path.
% Tests that run their own parfor over fibers (weighted regressions) are
% called serially per permutation and use the pool inside.
innerParallel = ismember(bundle.statfile, {'ea_explorer_stats_2sampleweightedlinreg', 'ea_explorer_stats_1sampleweightedlinreg'});
nWorkers = [];

% Fiber sets (coordinates once per run, not per analysis)
fibersets = rmfield(bundle.fibersets, 'valsin'); %#ok<NASGU>
save(fullfile(rundir, 'fibers.mat'), 'fibersets', '-v7.3');

settings = rmfield(bundle, 'fibersets');
shuffleInfo = {shuffles.info};
permsetFingerprint = permset.fingerprint;
permsetFile = opts.permsetFile;
checkReport = opts.checkReport;
runOptions = rmfield(opts, 'checkReport');
save(fullfile(rundir, 'run_info.mat'), 'settings', 'shuffleInfo', 'permsetFingerprint', ...
    'permsetFile', 'checkReport', 'runOptions', '-v7.3');

local_readme_start(rundir, bundle, permset, opts, shuffles);

% Small struct for ea_permtest_fibers_outcome, so workers do not receive
% the whole bundle (fiber coordinates)
lite = struct('mirrorsides', bundle.mirrorsides, 'Nsides', bundle.Nsides, ...
    'analyses', struct('gpatsel', {bundle.analyses.gpatsel}));

Nsides = bundle.Nsides;
statfile = bundle.statfile;
H0 = bundle.H0;
resultFiles = cell(1, numel(bundle.analyses));
seconds = zeros(1, numel(bundle.analyses));
compute = struct('label', {bundle.analyses.label}, 'path', '', 'reason', '', 'maxdiffR', NaN, 'maxdiffP', NaN);
usedNames = {};
for a = 1:numel(bundle.analyses)
    tA = tic;
    an = bundle.analyses(a);
    fs = bundle.fibersets(an.fiberset);
    Nfib = cellfun(@numel, fs.nonemptyidx);
    NpermA = size(shuffles(a).PermIdx, 2);

    fname = ['results_', regexprep(an.label, '[^\w-]', '_')];
    while ismember(fname, usedNames)
        fname = [fname, '_', num2str(a)]; %#ok<AGROW>
    end
    usedNames{end+1} = fname; %#ok<AGROW>
    resultFiles{a} = fullfile(rundir, [fname, '.mat']);

    fprintf('\n[%d/%d] %s: %d permutations, %s fibers\n', a, numel(bundle.analyses), an.label, NpermA, ...
        strjoin(arrayfun(@(s) sprintf('side %d: %d', s, Nfib(s)), 1:Nsides, 'Uni', 0), ', '));

    m = matfile(resultFiles{a}, 'Writable', true);
    m.label = an.label;
    m.outcomeCol = an.outcomeCol;
    m.outcome = bundle.outcomes(:, an.outcomeCol);
    m.patientselection = an.patientselection;
    m.gpatsel = an.gpatsel;
    m.fiberset = an.fiberset;
    m.Nperm = NpermA;
    m.PermIdx = shuffles(a).PermIdx;
    m.shuffleInfo = shuffles(a).info;
    m.statfile = statfile;

    % Real (unshuffled) statistic
    oc = ea_permtest_fibers_outcome(bundle.outcomes(:, an.outcomeCol), lite, a);
    for side = 1:Nsides
        R = nan(Nfib(side), 1);
        P = nan(Nfib(side), 1);
        if Nfib(side) > 0
            [R, P] = feval(statfile, fs.valsin{side}, oc{side}, H0);
        end
        m.(sprintf('fibIdx_side%d', side)) = fs.nonemptyidx{side};
        m.(sprintf('Remp_side%d', side)) = R(:);
        m.(sprintf('pemp_side%d', side)) = P(:);
        m.(sprintf('Rperm_side%d', side)) = nan(NpermA, Nfib(side), 'single');
        m.(sprintf('pperm_side%d', side)) = nan(NpermA, Nfib(side), 'single');
    end

    % Fast path for Spearman / Pearson, checked against the test function
    % on the real outcome and the first two shuffles
    Iperm = shuffles(a).Iperm;
    prep = cell(1, Nsides);
    useFast = opts.fastPath;
    reason = '';
    if ~opts.fastPath
        reason = 'switched off (opts.fastPath)';
    end
    dR = 0;
    dP = 0;
    for side = 1:Nsides
        if ~useFast || Nfib(side) == 0
            continue
        end
        prep{side} = ea_permtest_fastcorr_prepare(statfile, fs.valsin{side}, oc{side});
        if ~prep{side}.ok
            useFast = false;
            reason = prep{side}.reason;
            continue
        end
        nchk = min(2, NpermA);
        Ychk = oc{side};
        Rchk = m.(sprintf('Remp_side%d', side));
        Pchk = m.(sprintf('pemp_side%d', side));
        for j = 1:nchk
            ocj = ea_permtest_fibers_outcome(Iperm(:, j), lite, a);
            [r, q] = feval(statfile, fs.valsin{side}, ocj{side}, H0);
            Ychk = [Ychk, ocj{side}]; %#ok<AGROW>
            Rchk = [Rchk, r(:)]; %#ok<AGROW>
            Pchk = [Pchk, q(:)]; %#ok<AGROW>
        end
        [Rf, Pf] = ea_permtest_fastcorr(prep{side}, Ychk);
        dR = max([dR; abs(Rf(:) - Rchk(:))], [], 'omitnan');
        dP = max([dP; abs(Pf(:) - Pchk(:))], [], 'omitnan');
        if ~isequal(isnan(Rf), isnan(Rchk)) || ~isequal(isnan(Pf), isnan(Pchk)) || dR > 1e-10 || dP > 1e-10
            useFast = false;
            reason = sprintf('check against the test function failed on side %d (max |dR| %.1e, max |dP| %.1e)', side, dR, dP);
        end
    end
    if useFast
        compute(a).path = 'fast';
        compute(a).maxdiffR = dR;
        compute(a).maxdiffP = dP;
        fprintf('  Fast path (checked against %s: max |dR| %.1e, max |dP| %.1e)\n', statfile, dR, dP);
    else
        compute(a).path = 'test function';
        compute(a).reason = reason;
        fprintf('  Test function %s (fast path not used: %s)\n', statfile, reason);
        if isempty(nWorkers)
            nWorkers = ea_permtest_parpool(opts.maxWorkers);
        end
        if innerParallel
            nOuter = 0;
        else
            nOuter = nWorkers;
        end
    end

    % Shuffled statistics, in batches
    valsin = fs.valsin;
    maxNfib = max([Nfib, 1]);
    for b0 = 1:opts.batchSize:NpermA
        idx = b0:min(b0 + opts.batchSize - 1, NpermA);
        Iblk = Iperm(:, idx);
        if useFast
            Y = repmat({nan(numel(an.gpatsel), numel(idx))}, 1, Nsides);
            for j = 1:numel(idx)
                ocj = ea_permtest_fibers_outcome(Iblk(:, j), lite, a);
                for side = 1:Nsides
                    Y{side}(:, j) = ocj{side};
                end
            end
            for side = 1:Nsides
                if Nfib(side) > 0
                    [Rf, Pf] = ea_permtest_fastcorr(prep{side}, Y{side});
                    m.(sprintf('Rperm_side%d', side))(idx, 1:Nfib(side)) = single(Rf');
                    m.(sprintf('pperm_side%d', side))(idx, 1:Nfib(side)) = single(Pf');
                end
            end
        else
            Rb = nan(numel(idx), maxNfib, Nsides, 'single');
            Pb = nan(numel(idx), maxNfib, Nsides, 'single');
            parfor (j = 1:numel(idx), nOuter)
                ocj = ea_permtest_fibers_outcome(Iblk(:, j), lite, a); %#ok<PFBNS>
                rj = nan(1, maxNfib, Nsides, 'single');
                pj = nan(1, maxNfib, Nsides, 'single');
                for s = 1:Nsides
                    if Nfib(s) > 0
                        [r, q] = feval(statfile, valsin{s}, ocj{s}, H0); %#ok<PFBNS>
                        rj(1, 1:Nfib(s), s) = r(:)';
                        pj(1, 1:Nfib(s), s) = q(:)';
                    end
                end
                Rb(j, :, :) = rj;
                Pb(j, :, :) = pj;
            end
            for side = 1:Nsides
                if Nfib(side) > 0
                    m.(sprintf('Rperm_side%d', side))(idx, 1:Nfib(side)) = Rb(:, 1:Nfib(side), side);
                    m.(sprintf('pperm_side%d', side))(idx, 1:Nfib(side)) = Pb(:, 1:Nfib(side), side);
                end
            end
        end
        elapsed = toc(tA);
        fprintf('  %d/%d permutations (%.0f s, ~%.0f s left for this analysis)\n', idx(end), NpermA, ...
            elapsed, elapsed / idx(end) * (NpermA - idx(end)));
    end
    seconds(a) = toc(tA);
end

totalSeconds = toc(tStart);
if isempty(nWorkers)
    nWorkers = 0;
end
save(fullfile(rundir, 'run_info.mat'), 'seconds', 'totalSeconds', 'compute', 'nWorkers', '-append');
L = cell(1, numel(bundle.analyses));
for a = 1:numel(bundle.analyses)
    if strcmp(compute(a).path, 'fast')
        how = sprintf('fast path (checked: max |dR| %.1e, max |dP| %.1e)', compute(a).maxdiffR, compute(a).maxdiffP);
    elseif innerParallel
        how = sprintf('test function, %d worker(s) inside the test (%s)', max(nWorkers, 1), compute(a).reason);
    else
        how = sprintf('test function, %d worker(s) (%s)', max(nWorkers, 1), compute(a).reason);
    end
    L{a} = sprintf('- %s: %.0f s, %s -> `%s`', bundle.analyses(a).label, seconds(a), how, local_name(resultFiles{a}));
end
L{end+1} = sprintf('- Total: %.0f s (results do not depend on the computation path or the number of workers)', totalSeconds);
ea_permtest_readme_append(rundir, 'Permutations finished', L);

fprintf('\nDone in %.0f s. Results in:\n%s\n', totalSeconds, rundir);

run = struct('rundir', rundir, 'setupdir', setupdir, 'name', opts.name, 'resultFiles', {resultFiles}, ...
    'nWorkers', nWorkers, 'seconds', seconds, 'totalSeconds', totalSeconds, 'compute', compute);


function local_readme_start(rundir, bundle, permset, opts, shuffles)
L = {};
L{end+1} = sprintf('- Analysis: `%s`', opts.name);
L{end+1} = sprintf('- Tool: %s, explorer ID: `%s`, mode: %s', bundle.tool, bundle.ID, bundle.mode);
L{end+1} = sprintf('- Lead group folder: `%s` (%s)', bundle.groupdir, bundle.groupdirSource);
L{end+1} = sprintf('- Connectome: %s, fiber metric: %s', bundle.connectome, bundle.methodid);
L{end+1} = sprintf('- Statistic: %s (%s), H0: %s, test function: `%s`', bundle.stattest, bundle.statfamily, ...
    char(string(bundle.H0)), bundle.statfile);
L{end+1} = sprintf('- Coverage threshold: %g%% of patients, E-field threshold: %g', ...
    bundle.statsettings.connthreshold, bundle.statsettings.efieldthreshold);
L{end+1} = sprintf('- Mirrored E-fields: %d, selected patients: %d of %d', bundle.mirrorsides, ...
    numel(bundle.patientselection), bundle.Npatients);
if opts.stratifyByGroup
    L{end+1} = '- Shuffles: within cohort groups (stratified)';
else
    L{end+1} = '- Shuffles: across all selected patients (pooled)';
end
L{end+1} = sprintf('- Permutations: %d, permset fingerprint `%s`, file `%s`', permset.Nperm, ...
    permset.fingerprint, opts.permsetFile);
L{end+1} = '';
L{end+1} = '| Analysis | Patients | Missing | Fibers per side | Shuffle mode | Replaced | Jointness |';
L{end+1} = '|---|---|---|---|---|---|---|';
for a = 1:numel(bundle.analyses)
    an = bundle.analyses(a);
    fs = bundle.fibersets(an.fiberset);
    in = shuffles(a).info;
    k = shuffles(a).outcomePos;
    L{end+1} = sprintf('| %s | %d | %d | %s | %s | %d | %.2f |', an.label, numel(an.patientselection), ...
        in.NexcludedNaN(k), strjoin(cellfun(@(x) num2str(numel(x)), fs.nonemptyidx, 'Uni', 0), ' / '), ...
        in.mode, in.Nreplaced, in.jointness); %#ok<AGROW>
end
if isstruct(opts.checkReport) && isfield(opts.checkReport, 'ok')
    L{end+1} = '';
    if opts.checkReport.ok
        L{end+1} = '- Consistency check against the explorer: **passed**';
    else
        L{end+1} = '- Consistency check against the explorer: **FAILED**';
    end
    for a = 1:numel(opts.checkReport.analysis)
        sides = opts.checkReport.analysis(a).side;
        L{end+1} = sprintf('  - %s: compared %s fibers, max |diff| %s', opts.checkReport.analysis(a).label, ...
            strjoin(arrayfun(@(s) num2str(s.Ncompared), sides, 'Uni', 0), ' / '), ...
            strjoin(arrayfun(@(s) num2str(s.maxdiff), sides, 'Uni', 0), ' / ')); %#ok<AGROW>
    end
end
ea_permtest_readme_append(rundir, 'Fiber permutation run', L);


function s = local_name(f)
[~, n, e] = fileparts(f);
s = [n, e];


function opts = local_defaults(opts, defaults)
fn = fieldnames(defaults);
for i = 1:numel(fn)
    if ~isfield(opts, fn{i}) || (isempty(opts.(fn{i})) && ~isstruct(defaults.(fn{i})))
        opts.(fn{i}) = defaults.(fn{i});
    end
end

