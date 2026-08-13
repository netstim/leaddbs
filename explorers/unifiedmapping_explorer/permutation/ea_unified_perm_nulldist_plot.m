function [h, pval, exceedCount] = ea_unified_perm_nulldist_plot(nullvals, empval, xlabelstr, titlestr, tail, useTailApprox)
% Plots a permutation null distribution with the empirical value marked
% and its rank/p-value annotated in the title. Publication style: white
% background, no box, thin lines. Shared by the voxelwise-map-level plots
% (omnibus, max-statistic) and the out-of-sample prediction plot.
%
% tail: 'right' (nullvals >= empval, e.g. omnibus Q, positive max-stat)
%       'left'  (nullvals <= empval, e.g. negative max-stat)
%       'both'  (abs(nullvals) >= abs(empval), e.g. out-of-sample prediction R,
%                if the caller explicitly wants a two-sided test there)
% Defaults to 'right'.
%
% Callers in this feature deliberately do NOT collapse the voxelwise/
% max-stat/Eisenstein-Q tests to a single tail: sweet-spot (positive) and
% sour-spot (negative) directions are tracked and reported separately
% throughout, since a region where more field predicts WORSE outcomes is as
% scientifically informative as one where it predicts better. Only the
% out-of-sample prediction step (ea_unified_perm_sweetspot_predict) exposes
% `tail` as a caller-facing choice, defaulting to 'right'.
%
% useTailApprox (default false): if true and the rank-based p-value is
% already small, refine it with a Generalized Pareto Distribution tail fit
% (ea_unified_perm_palm_gpd_pval, ported from FSL PALM's palm_pareto.m) instead of
% leaving it capped at the permutation-count floor of 1/(Nperm+1). Opt-in
% and off by default so it never silently changes a p-value from a run
% that didn't ask for it. Falls back to the plain rank-based p-value if the
% tail doesn't fit a GPD well (checked internally), so this never
% fabricates a number. See ea_unified_perm_palm_gpd_pval.m for citations if used.
%
% NaN entries in nullvals are never counted as exceeding (comparisons
% against NaN are always false), but remain in the denominator -- a
% permutation that produced no defined statistic is evidence against the
% null, not missing data.

if ~exist('tail', 'var') || isempty(tail)
    tail = 'right';
end
if ~exist('useTailApprox', 'var') || isempty(useTailApprox)
    useTailApprox = false;
end

pal = ea_unified_perm_palette();

n = numel(nullvals);
switch tail
    case 'right'
        exceedCount = sum(nullvals >= empval);
    case 'left'
        exceedCount = sum(nullvals <= empval);
    case 'both'
        exceedCount = sum(abs(nullvals) >= abs(empval));
end
% The empirical value is itself one valid draw under the null (it's just
% one particular labeling, no more special than any permuted one if H0 is
% true), so it counts as a member of its own reference set: p = (b+1)/(m+1)
% rather than b/m. Without the +1, zero exceedances would report p = 0,
% which overstates certainty given only a finite number of permutations
% (Phipson & Smyth 2010, "Permutation p-values should never be zero").
pval = (exceedCount + 1) / (n + 1);

tailFitUsed = false;
if useTailApprox
    [pval, tailFitUsed] = ea_unified_perm_palm_gpd_pval(nullvals, empval, tail);
end

h = figure('Color', 'w');
histogram(nullvals, 'FaceColor', pal.blue, 'EdgeColor', 'none', 'FaceAlpha', 0.75);
hold on;
xline(empval, 'Color', pal.red, 'LineWidth', 2);
legend({'Null distribution', 'Empirical'}, 'Box', 'off', 'Location', 'best');
xlabel(xlabelstr);
ylabel('Number of permutations');
if tailFitUsed
    ptag = sprintf('p = %.4g (GPD tail fit)', pval);
else
    ptag = sprintf('p = %.4g', pval);
end
title({titlestr, sprintf('Empirical = %.3f, ranked %d of %d (%s)', empval, exceedCount, n, ptag)});
ea_unified_perm_style_axes(gca);
