function prep = ea_permtest_fastcorr_prepare(statfile, valsin, outcome)
% Fast path for the explorers' correlation tests (ea_explorer_stats_spearman,
% ea_explorer_stats_pearson) in permutation runs. These tests call
%   corr(valsin', outcome, 'Rows', 'pairwise', 'Type', <type>)
% once per permutation, which ranks every fiber again although only the
% outcome changes. Here the fiber values are ranked (Spearman) or
% standardised (Pearson) once; ea_permtest_fastcorr.m then computes R and p
% for many shuffled outcomes with one matrix product.
%
% Spearman reproduces corr bit for bit: with midranks (multiples of 0.5)
% the matrix product gives the sums of squared rank differences exactly,
% and R and p are then computed with corr's own formulas in corr's order.
% This matters because fibers with tied R values are common (discrete
% ranks, many zero E-fields) and later steps (display thresholds, spatial
% correlations of models) depend on these ties. corr's p-value rules:
% - n >= 10: t-approximation when there are ties, the AS 89 approximation
%   (Best & Roberts 1975) when there are none. Without missing scores corr
%   decides "ties" once for ALL fibers together (any tied fiber or tied
%   outcome); with missing scores it works fiber by fiber ('pairwise').
% - n < 10: corr enumerates all orderings (exact p) -- not reproduced, the
%   fast path is refused (the test function is used).
% Pearson (continuous values, ties practically impossible) agrees with corr
% up to rounding (~1e-15).
%
% statfile - the explorer's test function (bundle.statfile).
% valsin   - Nfibers x N fiber values, as passed to the test.
% outcome  - N x 1 real outcome, as passed to the test (NaN = missing).
%            Shuffled outcomes must have NaN in the same places (true for
%            the shuffles of ea_permtest_derive_shuffles.m).
%
% prep - struct: .ok (false: use the test function, .reason says why),
%        .type, .rows (patients with a score), .n, and
%        Spearman: .Xr (midranks), .sx2, .xadj, .X0, .dx, .tieBranch
%        Pearson:  .Zx (standardised values)

prep = struct('ok', false, 'reason', '', 'type', '', 'rows', [], 'n', 0);
switch statfile
    case 'ea_explorer_stats_spearman'
        prep.type = 'spearman';
    case 'ea_explorer_stats_pearson'
        prep.type = 'pearson';
    otherwise
        prep.reason = sprintf('no fast path for %s', statfile);
        return
end
if isempty(valsin)
    prep.reason = 'no fibers';
    return
end
if any(isnan(valsin(:)))
    prep.reason = 'missing fiber values';
    return
end

outcome = outcome(:);
rows = ~isnan(outcome);
n = sum(rows);
if strcmp(prep.type, 'spearman') && n < 10
    prep.reason = 'fewer than 10 patients (corr uses exact Spearman p-values)';
    return
end
if n < 3
    prep.reason = 'fewer than 3 patients';
    return
end

X = double(valsin(:, rows));
switch prep.type
    case 'spearman'
        [Xr, fiberTies, xadj] = ea_permtest_midranks(X);
        [~, outcomeTies] = ea_permtest_midranks(outcome(rows)');
        if any(~rows)
            % corr works fiber by fiber when scores are missing
            prep.tieBranch = fiberTies | outcomeTies;
        else
            prep.tieBranch = repmat(any(fiberTies) || outcomeTies, size(X, 1), 1);
        end
        prep.Xr = Xr;
        prep.sx2 = sum(Xr.^2, 2);
        prep.xadj = xadj;
        % centred ranks for corr's no-ties branch (Pearson of ranks)
        prep.X0 = Xr - (n + 1) / 2;
        prep.dx = vecnorm(prep.X0', 2, 1)';
    case 'pearson'
        X = X - mean(X, 2);
        prep.Zx = X ./ vecnorm(X, 2, 2); % constant fibers -> NaN, as in corr
end
prep.rows = rows;
prep.n = n;
prep.ok = true;
