function report = ea_permtest_check_disctract(obj, bundle, tol, compareAll)
% Verify that the permutation engine reproduces the fiber filtering
% explorer's own (unshuffled) fiber values: computes the statistic per
% fiber from the adapter's bundle, for every analysis (single tract, each
% subscore or each group), and compares it with ea_discfibers_calcstats(obj)
% on the fibers the explorer returns for that analysis. Run this once per
% analysis setup before trusting a permutation run.
%
% The explorer only returns the fibers its display settings show (e.g. the
% top 25%). With compareAll (default), the display settings are switched
% to "show every positive and negative fiber, no significance masking"
% for the duration of the check, so ALL analysed fibers are compared, not
% just the strongest ones.
%
% Everything the check changes on obj is restored afterwards (display
% settings, and obj.stats / obj.subscore, where ea_discfibers_calcstats
% stores fiber counts for the GUI), so obj is left as it was found.
%
% obj        - ea_disctract object.
% bundle     - from ea_permtest_adapter_disctract(obj). Built if omitted.
% tol        - maximum allowed absolute difference. Default 1e-6.
% compareAll - compare all analysed fibers (true, default) or only those
%              the explorer currently displays (false).
%
% report - struct with .ok and per analysis/side (report.analysis(a).side(s))
%          .Ncompared, .Nanalysed, .maxdiff, .nanmismatch, .outsideAnalysed,
%          .notReturned (analysed fibers with a finite, non-zero value that
%          the explorer did not return; compareAll only).

if ~exist('bundle', 'var') || isempty(bundle)
    bundle = ea_permtest_adapter_disctract(obj);
end
if ~exist('tol', 'var') || isempty(tol)
    tol = 1e-6;
end
if ~exist('compareAll', 'var') || isempty(compareAll)
    compareAll = true;
end

before = struct('stats', obj.stats, 'subscore', obj.subscore, ...
    'threshstrategy', obj.threshstrategy, 'showposamount', obj.showposamount, ...
    'shownegamount', obj.shownegamount, 'posvisible', obj.posvisible, ...
    'negvisible', obj.negvisible, 'showsignificantonly', obj.showsignificantonly);
restoreObj = onCleanup(@() local_restore(obj, before));
if compareAll
    % 'Fixed Amount' with an infinite amount keeps every positive and
    % every negative fiber (see ea_fibValThresh in ea_discfibers_calcstats.m)
    obj.threshstrategy = 'Fixed Amount';
    obj.showposamount = [inf inf];
    obj.shownegamount = [inf inf];
    obj.posvisible = 1;
    obj.negvisible = 1;
    obj.showsignificantonly = 0;
    if isfield(obj.subscore, 'special_case') && ~obj.subscore.special_case ...
            && isfield(obj.subscore, 'vis') && isfield(obj.subscore, 'posvisible')
        % split modes use per-subscore/per-group display settings
        obj.subscore.posvisible(:) = 1;
        obj.subscore.negvisible(:) = 1;
        obj.subscore.vis.showposamount(:) = inf;
        obj.subscore.vis.shownegamount(:) = inf;
    end
end
[vals, ~, usedidx] = ea_discfibers_calcstats(obj);
clear restoreObj

report = struct('ok', true, 'analysis', struct([]));
for a = 1:numel(bundle.analyses)
    an = bundle.analyses(a);
    fs = bundle.fibersets(an.fiberset);
    outcomein = ea_permtest_fibers_outcome(bundle.outcomes(:, an.outcomeCol), bundle, a);
    report.analysis(a).label = an.label;
    for side = 1:bundle.Nsides
        Remp = nan(fs.Nfibers{side}, 1);
        if ~isempty(fs.valsin{side})
            valsout = feval(bundle.statfile, fs.valsin{side}, outcomein{side}, bundle.H0);
            Remp(fs.nonemptyidx{side}) = valsout;
        end

        if size(vals, 1) >= an.key && ~isempty(vals{an.key, side})
            live = vals{an.key, side}(:);
            u = usedidx{an.key, side}(:);
        else
            live = zeros(0, 1);
            u = zeros(0, 1);
        end
        outsideAnalysed = sum(~ismember(u, fs.nonemptyidx{side}));
        nanmismatch = sum(isnan(Remp(u)) ~= isnan(live));
        d = abs(live - Remp(u));
        maxdiff = max([0; d(~isnan(d))]);
        if compareAll
            % every fiber with a finite, non-zero value should be returned
            notReturned = sum(isfinite(Remp) & Remp ~= 0 & ~ismember((1:numel(Remp))', u));
        else
            notReturned = 0;
        end

        r = struct('Ncompared', numel(u), 'Nanalysed', numel(fs.nonemptyidx{side}), ...
            'maxdiff', maxdiff, 'nanmismatch', nanmismatch, 'outsideAnalysed', outsideAnalysed, ...
            'notReturned', notReturned);
        r.ok = maxdiff <= tol && nanmismatch == 0 && outsideAnalysed == 0 && notReturned == 0;
        report.analysis(a).side(side) = r;
        report.ok = report.ok && r.ok;

        fprintf(['%s, side %d: %d analysed fibers, %d compared with explorer, max |diff| = %g, ', ...
            'NaN mismatches = %d, outside analysed set = %d, not returned by explorer = %d\n'], ...
            an.label, side, r.Nanalysed, r.Ncompared, maxdiff, nanmismatch, outsideAnalysed, notReturned);
        if r.Ncompared == 0
            ea_cprintf('CmdWinWarnings', '  (nothing to compare: the explorer currently shows no fibers for %s, side %d)\n', an.label, side);
        end
    end
end

if report.ok
    fprintf('Consistency check passed: permutation engine reproduces the explorer''s fiber values.\n');
else
    ea_cprintf('CmdWinWarnings', ['Consistency check FAILED: the permutation engine does not reproduce ', ...
        'the explorer''s fiber values -- do not use permutation results for this analysis until resolved.\n']);
end


function local_restore(obj, before)
fn = fieldnames(before);
for i = 1:numel(fn)
    obj.(fn{i}) = before.(fn{i});
end
