function [pval, exceedCount, tailFitUsed] = ea_permtest_nulldist_panel(ax, nullvals, empval, xlabelstr, titlestr, useTailApprox, thresh, threshLabel)
% Draw one permutation null distribution into axes ax, with the real
% (empirical) value marked and its rank-based p-value in the title.
% Publication style as in the sweetspot permutation plots (white
% background, no box, thin lines; see ea_permtest_style_axes).
% Right tail: large values are extreme (omnibus Q, positive max-statistic).
%
% ax            - axes to draw into.
% nullvals      - Nperm x 1 null distribution.
% empval        - real value (e.g. Q of the real map, or the highest real R).
% xlabelstr, titlestr - labels.
% useTailApprox - refine small p-values with the GPD tail fit
%                 (ea_permtest_gpd_pval). Default false.
% thresh        - optional value to mark with a dashed line (e.g. the FWE
%                 threshold of a max-statistic null), with threshLabel.
%
% p = (count + 1) / (Nperm + 1), count = permutations >= empval (NaN never
% counts as exceeding but stays in the denominator), as in the sweetspot
% permutation code.

if ~exist('useTailApprox', 'var') || isempty(useTailApprox)
    useTailApprox = false;
end
pal = ea_permtest_palette();

nullvals = nullvals(:);
n = numel(nullvals);
exceedCount = sum(nullvals >= empval);
pval = (exceedCount + 1) / (n + 1);
tailFitUsed = false;
if useTailApprox
    [pval, tailFitUsed] = ea_permtest_gpd_pval(nullvals, empval, 'right');
end

histogram(ax, nullvals(isfinite(nullvals)), 'FaceColor', pal.blue, 'EdgeColor', 'none', 'FaceAlpha', 0.75);
hold(ax, 'on');
hEmp = xline(ax, empval, 'Color', pal.red, 'LineWidth', 2);
handles = hEmp;
labels = {'Real data'};
if exist('thresh', 'var') && ~isempty(thresh) && isfinite(thresh)
    hThr = xline(ax, thresh, '--', 'Color', pal.ink, 'LineWidth', 1);
    handles(end+1) = hThr;
    labels{end+1} = threshLabel;
end
legend(ax, handles, labels, 'Box', 'off', 'Location', 'best');
xlabel(ax, xlabelstr);
ylabel(ax, 'Permutations');
if tailFitUsed
    ptag = sprintf('p = %.4g (GPD tail fit)', pval);
else
    ptag = sprintf('p = %.4g', pval);
end
title(ax, {titlestr, sprintf('real = %.3g, exceeded by %d of %d (%s)', empval, exceedCount, n, ptag)}, ...
    'FontWeight', 'normal', 'Interpreter', 'none');
ea_permtest_style_axes(ax);
