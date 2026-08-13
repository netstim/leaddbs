function ea_unified_perm_sweetspot_predict_plot(predpermtestFile, responselabel)
% Regenerates the two prediction-validation figures (empirical correlation
% scatter + out-of-sample null distribution) from a saved results.mat file
% (produced by ea_unified_perm_sweetspot_predict), without rerunning the
% prediction. Standalone -- only needs the file path, no live
% ea_unifiedmapping object. Called automatically at the end of
% ea_unified_perm_sweetspot_predict, and safe to call again later to
% regenerate the figures from a saved file.
%
% responselabel (optional): x-axis label for the correlation scatter plot.
% Defaults to the label saved at prediction time, or 'Response Variable' if
% that's unavailable.
%
% Both figures, plus a README section describing them, are saved into the
% same folder predpermtestFile lives in.

pred = load(predpermtestFile);
outdir = [fileparts(predpermtestFile), filesep];
[~, basefname] = fileparts(fileparts(predpermtestFile)); % predpermtestFile lives in its own run folder -- use that folder's name to tag exports

if ~exist('responselabel', 'var') || isempty(responselabel)
    if isfield(pred, 'responsevarlabel') && ~isempty(pred.responsevarlabel)
        responselabel = pred.responsevarlabel;
    else
        responselabel = 'Response Variable';
    end
end

if isfield(pred, 'tail') && ~isempty(pred.tail)
    tail = pred.tail;
else
    tail = 'both'; % older/foreign files that predate the tail field
end

% Figure 1: empirical Ihat vs. the test cohort's real, unpermuted responsevar.
hScatterSaved = false;
if isfield(pred, 'Iemp') && isfield(pred, 'Ihatemp')
    h1 = ea_corrbox(pred.Iemp, pred.Ihatemp, pred.pperm, ...
        {sprintf('Out-of-Sample Prediction (%s)', pred.testID), responselabel, 'Predicted Score (Ihat)'});
    try % ea_corrbox is a shared, gramm-based utility -- best-effort cleanup, non-fatal if its internal structure doesn't cooperate
        set(h1, 'Color', 'w');
        axesInH1 = findall(h1, 'Type', 'axes');
        for a = 1:numel(axesInH1)
            set(axesInH1(a), 'Color', 'w');
            box(axesInH1(a), 'off');
        end
    end
    saveas(h1, fullfile(outdir, sprintf('%s_scatter.png', basefname)));
    hScatterSaved = true;
else
    warning('This results.mat file has no Iemp/Ihatemp -- skipping the correlation scatter (Figure 1). Showing the null-distribution plot (Figure 2) only.');
end

% Figure 2: null distribution of out-of-sample R, empirical R marked, rank annotated.
h2 = ea_unified_perm_nulldist_plot(pred.Rpredperm, pred.Rpredemp, ...
    sprintf('Out-of-sample prediction R (%s)', pred.corrType), ...
    sprintf('Out-of-Sample Permutation Test (%d/%d permutations produced no defined prediction)', pred.nNaNperm, numel(pred.Rpredperm)), ...
    tail);
saveas(h2, fullfile(outdir, sprintf('%s_nulldist.png', basefname)));

readmeBody = sprintf([ ...
    'Out-of-sample prediction validation for %s, against training permtest %s.\n\n', ...
    'Settings used for this run:\n', ...
    '  sigMode  = %s\n', ...
    '  tail     = %s\n', ...
    '  Rpredemp = %.4f (parametric p = %.4g)\n', ...
    '  pperm    = %.4g (%d/%d permutations as extreme as Rpredemp, tail=''%s'')\n\n', ...
    'Files added by this step:\n\n', ...
    '  %s_scatter.png%s\n', ...
    '    This test cohort''s real, unpermuted response scores vs. the score predicted\n', ...
    '    by projecting the TRAINING cohort''s empirical (unpermuted) map onto this\n', ...
    '    cohort''s own patients. Rpredemp above is the correlation shown here.\n\n', ...
    '  %s_nulldist.png\n', ...
    '    Null distribution of out-of-sample prediction R: for each of the training\n', ...
    '    permtest''s Nperm permuted maps, the same projection+prediction was repeated\n', ...
    '    using this cohort''s real scores (never permuted -- only the training map\n', ...
    '    changes per draw). Answers "does the trained map predict THIS cohort''s\n', ...
    '    outcomes better than a map built from randomly-labeled training data would."\n', ...
    '    pperm above is read off this distribution.\n'], ...
    pred.testID, pred.trainPermtestFile, pred.sigMode, tail, pred.Rpredemp, pred.Rpredemp_pval, ...
    pred.pperm, pred.exceedCount, numel(pred.Rpredperm), tail, ...
    basefname, ea_unified_perm_ternary(hScatterSaved, '', ' (not generated -- Iemp/Ihatemp unavailable)'), basefname);
ea_unified_perm_readme_append(outdir, 'Out-of-sample prediction (ea_unified_perm_sweetspot_predict)', readmeBody);
end

function out = ea_unified_perm_ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end
