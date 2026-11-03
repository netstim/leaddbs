function report = ea_permtest_verify_disctract(obj, rundir, perms, tol)
% Independent check of a finished fiber permutation run against the fiber
% filtering explorer itself.
%
% Part 1 - shuffles: the saved PermIdx of every analysis is re-derived from
%   the permset file the run used and must be identical; every shuffle may
%   only rearrange the analysis' own selected patients' scores, and
%   patients missing that score must keep it (NaN stays in place).
% Part 2 - permuted maps: for the chosen permutations, the shuffled scores
%   are put into the explorer and ea_discfibers_calcstats (the explorer's
%   own code) recomputes the fiber maps; these are compared fiber by fiber
%   with the saved Rperm rows of the run.
%   - single tract / group mode: the shuffled outcome is passed to
%     ea_discfibers_calcstats through its permutation argument;
%   - subscore mode: ea_discfibers_calcstats ignores that argument and
%     reads obj.subscore.vars, so the subscores are temporarily replaced
%     by their shuffled versions.
%   Display settings are temporarily set to show all fibers (as in the
%   consistency check). Everything changed on obj is restored afterwards.
%
% obj    - the ea_disctract object the run was computed from (same
%          settings; checked).
% rundir - the run folder.
% perms  - permutation numbers to recompute. Default: first, two random,
%          last.
% tol    - maximum allowed absolute difference in R. Default 1e-5 (the
%          run stores shuffled R in single precision, ~7 significant digits).
%
% report - struct: .ok, .shuffles(a) and .maps(k) with per permutation,
%          analysis and side: Ncompared, maxdiff, Nmismatch, NnotReturned.

if ~exist('tol', 'var') || isempty(tol)
    tol = 1e-5;
end
ri = load(fullfile(rundir, 'run_info.mat'));
settings = ri.settings;
mapset = ea_permtest_read_fiberrun(rundir);
Nperm = mapset.maps(1).Nperm;
if ~exist('perms', 'var') || isempty(perms)
    rs = RandStream('mt19937ar', 'Seed', 1);
    perms = unique([1, sort(randperm(rs, Nperm, 2)), Nperm]);
end

%% The explorer must be in the same state as when the run was computed
same = strcmp(obj.multitractmode, settings.mode) && isequal(obj.statsettings, settings.statsettings) ...
    && isequal(obj.patientselection(:)', settings.patientselection) && logical(obj.mirrorsides) == settings.mirrorsides ...
    && strcmp(obj.connectome, settings.connectome);
if ~same
    ea_error(['The explorer settings differ from the run (mode, statistics settings, patient selection, ', ...
        'mirroring or connectome) -- load the explorer state the run was computed from.']);
end

report = struct('ok', true, 'perms', perms, 'shuffles', struct([]), 'maps', struct([]));

%% Part 1: shuffles
fprintf('\nPart 1: shuffles\n');
P = load(ri.permsetFile, 'permset');
if ~strcmp(P.permset.fingerprint, ri.permsetFingerprint)
    ea_error('The permset file does not match the fingerprint recorded in the run.');
end
bundle = ea_permtest_adapter_disctract(obj, struct('groupdir', settings.groupdir));
stratify = isfield(ri.runOptions, 'stratifyByGroup') && ri.runOptions.stratifyByGroup;
reshuf = ea_permtest_bundle_shuffles(P.permset, bundle, stratify);
for a = 1:numel(settings.analyses)
    m = matfile(local_file(rundir, settings.analyses(a).label));
    PermIdx = m.PermIdx;
    I = m.outcome;
    sel = settings.analyses(a).patientselection;
    own = sel(~isnan(I(sel)));
    others = setdiff(1:numel(I), own);
    identical = isequal(PermIdx, reshuf(a).PermIdx);
    onlyOwn = all(ismember(PermIdx(own, :), own), 'all');
    fixed = all(PermIdx(others, :) == others(:), 'all');
    isPerm = all(sort(PermIdx(own, :), 1) == sort(own(:)), 'all');
    r = struct('label', settings.analyses(a).label, 'identicalToPermset', identical, ...
        'onlyOwnPatients', onlyOwn, 'othersFixed', fixed, 'eachShuffleIsPermutation', isPerm);
    r.ok = identical && onlyOwn && fixed && isPerm;
    report.shuffles = [report.shuffles, r];
    report.ok = report.ok && r.ok;
    fprintf('  %-32s re-derived identical: %d | only own patients: %d | missing/unselected fixed: %d | each column a permutation: %d\n', ...
        r.label, identical, onlyOwn, fixed, isPerm);
end

%% Part 2: permuted maps recomputed by the explorer
fprintf('\nPart 2: permuted maps recomputed with ea_discfibers_calcstats, permutations %s\n', mat2str(perms));
before = struct('stats', obj.stats, 'subscore', obj.subscore, 'responsevar', obj.responsevar, ...
    'threshstrategy', obj.threshstrategy, 'showposamount', obj.showposamount, ...
    'shownegamount', obj.shownegamount, 'posvisible', obj.posvisible, ...
    'negvisible', obj.negvisible, 'showsignificantonly', obj.showsignificantonly);
restoreObj = onCleanup(@() local_restore(obj, before));
obj.threshstrategy = 'Fixed Amount';
obj.showposamount = [inf inf];
obj.shownegamount = [inf inf];
obj.posvisible = 1;
obj.negvisible = 1;
obj.showsignificantonly = 0;
if isfield(obj.subscore, 'special_case') && ~obj.subscore.special_case ...
        && isfield(obj.subscore, 'vis') && isfield(obj.subscore, 'posvisible')
    obj.subscore.posvisible(:) = 1;
    obj.subscore.negvisible(:) = 1;
    obj.subscore.vis.showposamount(:) = inf;
    obj.subscore.vis.shownegamount(:) = inf;
end

subscoreMode = strcmp(settings.mode, 'Split & Color By Subscore');
for p = perms
    if subscoreMode
        % all subscores replaced by their shuffled versions at once
        for a = 1:numel(settings.analyses)
            m = matfile(local_file(rundir, settings.analyses(a).label));
            I = m.outcome;
            PermIdx = m.PermIdx(:, p);
            obj.subscore.vars{settings.analyses(a).key} = I(PermIdx);
        end
        [vals, ~, usedidx] = ea_discfibers_calcstats(obj);
        for a = 1:numel(settings.analyses)
            report = local_compare(report, rundir, settings, mapset, a, p, vals, usedidx, tol);
        end
    else
        for a = 1:numel(settings.analyses)
            m = matfile(local_file(rundir, settings.analyses(a).label));
            I = m.outcome;
            PermIdx = m.PermIdx(:, p);
            [vals, ~, usedidx] = ea_discfibers_calcstats(obj, obj.patientselection, I(PermIdx));
            report = local_compare(report, rundir, settings, mapset, a, p, vals, usedidx, tol);
        end
    end
end
clear restoreObj

if report.ok
    fprintf('\nVerification passed: shuffles re-derive identically, and the explorer reproduces the saved permuted maps.\n');
else
    ea_cprintf('CmdWinWarnings', '\nVerification FAILED -- see the lines above.\n');
end


function report = local_compare(report, rundir, settings, mapset, a, p, vals, usedidx, tol)
an = settings.analyses(a);
m = matfile(local_file(rundir, an.label));
for side = 1:settings.Nsides
    fibIdx = m.(sprintf('fibIdx_side%d', side));
    if isempty(fibIdx)
        continue
    end
    stored = double(m.(sprintf('Rperm_side%d', side))(p, :))';
    if size(vals, 1) >= an.key && ~isempty(vals{an.key, side})
        live = vals{an.key, side}(:);
        u = usedidx{an.key, side}(:);
    else
        live = zeros(0, 1);
        u = zeros(0, 1);
    end
    [inRun, pos] = ismember(u, fibIdx);
    d = abs(live(inRun) - stored(pos(inRun)));
    returned = false(size(fibIdx));
    returned(pos(inRun)) = true;
    notReturned = sum(isfinite(stored) & stored ~= 0 & ~returned);
    r = struct('perm', p, 'label', an.label, 'side', side, 'Ncompared', nnz(inRun), ...
        'maxdiff', max([0; d]), 'Nmismatch', sum(d > tol), 'NoutsideRun', sum(~inRun), 'NnotReturned', notReturned);
    r.ok = r.Nmismatch == 0 && r.NoutsideRun == 0 && notReturned == 0 && r.Ncompared > 0;
    report.maps = [report.maps, r];
    report.ok = report.ok && r.ok;
    fprintf('  perm %4d  %-32s side %d: %5d fibers compared, max |diff| = %.2g, mismatches = %d, not in run = %d, not returned = %d\n', ...
        p, an.label, side, r.Ncompared, r.maxdiff, r.Nmismatch, r.NoutsideRun, notReturned);
end


function f = local_file(rundir, label)
f = fullfile(rundir, ['results_', regexprep(label, '[^\w-]', '_'), '.mat']);
if ~isfile(f)
    files = dir(fullfile(rundir, 'results_*.mat'));
    for i = 1:numel(files)
        if strcmp(matfile(fullfile(rundir, files(i).name)).label, label)
            f = fullfile(rundir, files(i).name);
            return
        end
    end
    ea_error(['No results file for ', label]);
end


function local_restore(obj, before)
fn = fieldnames(before);
for i = 1:numel(fn)
    obj.(fn{i}) = before.(fn{i});
end
