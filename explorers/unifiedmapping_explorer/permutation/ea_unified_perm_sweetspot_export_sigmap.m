function fnames = ea_unified_perm_sweetspot_export_sigmap(permtestFile, alphalevel, posvisible, negvisible, outdir)
% Standalone, manually-invoked utility: exports the empirical R-map from a
% saved permtest results.mat file, masked to voxels significant at
% pemp <= alphalevel -- this reproduces exactly the significance mask
% ea_unified_perm_sweetspot_predict applies to the empirical training map
% before out-of-sample prediction (see its Remp_train/nonsig logic), plus
% the positive/negative visibility mask applied at prediction time (see
% ea_unified_perm_sweetspot_predictfromvals: posvisible=0 blanks R>0,
% negvisible=0 blanks R<0).
%
% Not called automatically by anything else in this feature -- a tool for
% you to inspect/export a specific significance mask by hand.
%
% permtestFile - path to a results.mat saved by ea_unified_perm_sweetspot.
% alphalevel   - uncorrected significance cutoff on pemp (default 0.05).
% posvisible   - include positive ("sweet spot") voxels (default 1).
% negvisible   - include negative ("sour spot") voxels (default 1).
% outdir       - output folder (default: same folder as permtestFile, i.e.
%                that run's own folder).
%
% Writes <run>_empirical_sig_p<alpha>_r/_l.nii and returns the written file
% paths. Appends a README section to that folder describing what was added.

if ~exist('alphalevel', 'var') || isempty(alphalevel)
    alphalevel = 0.05;
end
if ~exist('posvisible', 'var') || isempty(posvisible)
    posvisible = 1;
end
if ~exist('negvisible', 'var') || isempty(negvisible)
    negvisible = 1;
end
if ~exist('outdir', 'var') || isempty(outdir)
    outdir = [fileparts(permtestFile), filesep];
end

S = load(permtestFile, 'Remp', 'pemp', 'space');
nsides = numel(S.space);
[~, basefname] = fileparts(fileparts(permtestFile)); % permtestFile lives in its own run folder -- use that folder's name to tag exports

sigMap = cell(1, nsides);
for side = 1:nsides
    Remp = S.Remp{side};
    pemp = S.pemp{side};

    nonsig = isnan(pemp) | pemp > alphalevel;
    v = Remp;
    v(nonsig) = nan;

    if ~posvisible
        v(v > 0) = nan;
    end
    if ~negvisible
        v(v < 0) = nan;
    end

    sigMap{side} = v;
    fprintf('Side %d: %d/%d voxels significant at p<=%.3g (posvisible=%d, negvisible=%d).\n', ...
        side, sum(~isnan(v)), numel(v), alphalevel, posvisible, negvisible);
end

fnames = ea_unified_perm_vals2nii(S.space, sigMap, outdir, sprintf('%s_empirical_sig_p%.3g', basefname, alphalevel));

readmeBody = sprintf([ ...
    'Manually-exported empirical significance mask for %s.\n\n', ...
    'Settings used:\n', ...
    '  alphalevel = %.3g (uncorrected, on pemp -- the empirical map''s own parametric p)\n', ...
    '  posvisible = %d\n', ...
    '  negvisible = %d\n\n', ...
    'Files added by this step:\n\n', ...
    '  %s_empirical_sig_p%.3g_r/_l.nii\n', ...
    '    The empirical (real, unpermuted) R-map, masked to voxels with pemp<=alphalevel\n', ...
    '    (this is the SAME uncorrected, per-voxel test as *_voxelwise_uncorrected_*\n', ...
    '    from ea_unified_perm_sweetspot_stats -- exported here on demand, independent\n', ...
    '    of whether that step has been run for this folder).\n'], ...
    permtestFile, alphalevel, posvisible, negvisible, basefname, alphalevel);
ea_unified_perm_readme_append(outdir, 'Manual empirical significance export (ea_unified_perm_sweetspot_export_sigmap)', readmeBody);
end
