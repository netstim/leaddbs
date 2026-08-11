function results = ea_unified_perm_sweetspot_tfce(permtestFile, E, H, conn, deltah, NpermTFCE, maxWorkers, alpha, chunkSize)
% CITE IF USED: applies PALM (Permutation Analysis of Linear Models)'s TFCE
% method (Threshold-Free Cluster Enhancement) to a saved sweetspot
% permutation-test results file, via ea_unified_perm_palm_tfce.m (a port of
% PALM's palm_tfce.m). If you report results from this function, cite:
%   Winkler AM, Ridgway GR, Douaud G, Nichols TE, Smith SM (2016). Faster
%   permutation inference in brain imaging. NeuroImage 141:502-516. (PALM)
%   Smith SM, Nichols TE (2009). Threshold-free cluster enhancement:
%   addressing problems of smoothing, threshold dependence and
%   localisation in cluster inference. NeuroImage 44:83-98. (TFCE itself)
%
% Reads a results.mat file (produced by ea_unified_perm_sweetspot) and
% computes a TFCE-enhanced, max-based FWE-corrected significance map --
% a fourth, complementary analysis alongside the three already computed by
% ea_unified_perm_sweetspot_stats.m (uncorrected voxelwise, Eisenstein
% omnibus, max-statistic). Positive ("sweet spot") and negative ("sour
% spot") directions are enhanced and thresholded separately, same
% convention as the rest of this feature: TFCE is one-sided by
% construction (see ea_unified_perm_palm_tfce.m), so the negative direction is
% obtained by enhancing the NEGATED map, not by any special-cased logic
% here.
%
% Correction logic (standard PALM/randomise usage of TFCE, not just the
% enhancement step alone): for each permutation, take the single largest
% TFCE value anywhere in the map -- this builds a null distribution
% structurally identical to ea_unified_perm_sweetspot_stats.m's max-stat
% test (same accumulate-a-running-max idea, same (count+1)/(Nperm+1)
% p-value convention), just applied to TFCE-enhanced values instead of raw
% R. Every voxel's empirical TFCE value is then ranked against this SAME
% whole-map-max null distribution -- this is what makes the result
% multiple-comparison corrected, not the enhancement step by itself.
%
% permtestFile - path to a results.mat saved by ea_unified_perm_sweetspot.
% E, H, conn, deltah - TFCE parameters, passed straight through to
%              ea_unified_perm_palm_tfce.m. Defaults 0.5 / 2 / 6 / 0(='auto'),
%              matching PALM's own defaults.
% NpermTFCE  - how many of the file's saved permutations to use for TFCE.
%              Defaults to ALL of them (the file's own Nperm) -- TFCE here
%              reuses the identical bwconncomp call PALM itself uses, so a
%              full-scale run should be feasible in a broadly comparable
%              timeframe to a native PALM run on similar data (PALM's own
%              TFCE loop is plain serial MATLAB with no parfor/parpool
%              anywhere in its source -- confirmed by inspection -- so the
%              parfor below should in principle make this faster than
%              PALM's own single-threaded runtime on comparable hardware).
%              Lower it if you want a faster, lower-resolution first look.
% maxWorkers - cap on parpool workers (ea_unified_perm_cap_parpool.m).
%              Default 3. Each worker here only needs one Nvoxels-sized row
%              (a few MB) plus a handful of scalars -- much lower memory
%              risk per worker than ea_unified_perm_sweetspot.m's own
%              parfor, which has to hold a full copy of the model object.
% alpha      - significance threshold for the exported thresholded map.
%              Default 0.05.
% chunkSize  - voxels processed per chunk when computing the final
%              per-voxel p-values against the (small, NpermTFCE-long) null
%              distribution -- same memory-safety technique as
%              ea_unified_perm_sweetspot_stats.m, smaller default (5000)
%              since this step's per-chunk temporary is NpermTFCE x
%              chunkSize for two directions at once. Default 5000.
%
% Outputs, all written into the same folder permtestFile lives in:
%   *_tfce_p<alpha>_r/_l.nii - thresholded significance map (FWE-corrected).
%   *_tfce_pos.png / *_tfce_neg.png - null-distribution plots (max TFCE
%     per permutation, empirical whole-map max TFCE marked).
%   README.txt (appended) - settings used, citation, explicit statement
%     that this map IS multiple-comparison-corrected.
%
% Returns a struct with everything computed ({side}-indexed).

if ~exist('E','var') || isempty(E)
    E = 0.5;
end
if ~exist('H','var') || isempty(H)
    H = 2;
end
if ~exist('conn','var') || isempty(conn)
    conn = 6;
end
if ~exist('deltah','var') || isempty(deltah)
    deltah = 0;
end
if ~exist('maxWorkers','var') || isempty(maxWorkers)
    maxWorkers = 3;
end
if ~exist('alpha','var') || isempty(alpha)
    alpha = 0.05;
end
if ~exist('chunkSize','var') || isempty(chunkSize)
    chunkSize = 5000;
end

meta = load(permtestFile, 'space', 'Remp', 'Nperm', 'nsides');
mFile = matfile(permtestFile);
outdir = [fileparts(permtestFile), filesep];
[~, basefname] = fileparts(fileparts(permtestFile)); % permtestFile lives in its own run folder -- use that folder's name to tag exports

nsides = meta.nsides;

if ~exist('NpermTFCE','var') || isempty(NpermTFCE)
    NpermTFCE = meta.Nperm;
end
if NpermTFCE > meta.Nperm
    ea_error(sprintf('NpermTFCE (%d) cannot exceed the %d permutations actually saved in %s.', NpermTFCE, meta.Nperm, permtestFile));
end

results = struct();
results.sourcefile = permtestFile;
results.E = E; results.H = H; results.conn = conn; results.deltah = deltah;
results.NpermTFCE = NpermTFCE;
results.alpha = alpha;

tfceThreshMap = cell(1, nsides);
hTfcePos = gobjects(1, nsides);
hTfceNeg = gobjects(1, nsides);

for side = 1:nsides
    dim = meta.space{side}.dim(1:3);
    Remp = meta.Remp{side};
    Nvox = numel(Remp);
    tag = sprintf('side%d', side);

    fprintf('%s: computing empirical TFCE (positive and negative directions)...\n', tag);
    tfceEmpPos = ea_unified_perm_palm_tfce(Remp, 'volume', dim, E, H, conn, deltah);
    tfceEmpNeg = ea_unified_perm_palm_tfce(-Remp, 'volume', dim, E, H, conn, deltah); % TFCE is one-sided by construction -- negative ("sour spot") direction is enhancing the negated map, not a special-cased branch

    RvarName = sprintf('Rperm_side%d', side);

    maxPosPerm = nan(NpermTFCE, 1);
    maxNegPerm = nan(NpermTFCE, 1);

    ea_unified_perm_cap_parpool(maxWorkers);

    fprintf('%s: running TFCE on %d permutations...\n', tag, NpermTFCE);
    parfor p = 1:NpermTFCE
        row = double(mFile.(RvarName)(p,:))'; %#ok<PFBNS> -- mFile/RvarName are broadcast, read-only inside parfor
        tp = ea_unified_perm_palm_tfce(row, 'volume', dim, E, H, conn, deltah);
        tn = ea_unified_perm_palm_tfce(-row, 'volume', dim, E, H, conn, deltah);
        maxPosPerm(p) = max(tp, [], 'omitnan');
        maxNegPerm(p) = max(tn, [], 'omitnan');
    end

    % Per-voxel FWE-corrected p-value: every voxel is ranked against the
    % SAME whole-map-max null distribution (not its own, unlike the
    % uncorrected voxelwise test in _stats.m) -- this is what makes TFCE's
    % result multiple-comparison corrected. Chunked over voxels so the
    % NpermTFCE x chunkSize broadcast temporary stays bounded regardless
    % of Nvoxels/NpermTFCE.
    pTfcePos = nan(1, Nvox);
    pTfceNeg = nan(1, Nvox);
    for vStart = 1:chunkSize:Nvox
        vEnd = min(vStart+chunkSize-1, Nvox);
        vIdx = vStart:vEnd;
        countPos = sum(maxPosPerm >= tfceEmpPos(vIdx)', 1); % transpose: indexing a column vector (tfceEmpPos) keeps its column orientation regardless of vIdx's shape, but maxPosPerm >= ... needs a row here to broadcast into an NpermTFCE x chunkSize comparison
        countNeg = sum(maxNegPerm >= tfceEmpNeg(vIdx)', 1);
        pTfcePos(vIdx) = (countPos + 1) / (NpermTFCE + 1);
        pTfceNeg(vIdx) = (countNeg + 1) / (NpermTFCE + 1);
    end
    pTfcePos(isnan(tfceEmpPos)) = nan;
    pTfceNeg(isnan(tfceEmpNeg)) = nan;

    survivePos = pTfcePos <= alpha;
    surviveNeg = pTfceNeg <= alpha;
    thisThreshMap = nan(size(Remp));
    thisThreshMap(survivePos') = Remp(survivePos');
    thisThreshMap(surviveNeg') = Remp(surviveNeg');
    tfceThreshMap{side} = thisThreshMap;

    results.tfceEmpPos{side} = tfceEmpPos;
    results.tfceEmpNeg{side} = tfceEmpNeg;
    results.pTfcePos{side} = pTfcePos;
    results.pTfceNeg{side} = pTfceNeg;
    results.maxPosPerm{side} = maxPosPerm;
    results.maxNegPerm{side} = maxNegPerm;
    results.thresholdedMap{side} = thisThreshMap;

    % TFCE-enhanced values are non-negative by construction regardless of
    % original R sign (extent^E * h^H, both directions enhanced off their
    % own non-negative h sweep) -- both null distributions use tail
    % 'right' (larger = more extreme), unlike ea_unified_perm_sweetspot_stats.m's
    % maxstat test where the raw R-based negative tail is genuinely
    % negative-valued and needs tail 'left'.
    [hTfcePos(side), results.pOmnibusPos{side}] = ea_unified_perm_nulldist_plot( ...
        maxPosPerm, max(tfceEmpPos, [], 'omitnan'), 'Max TFCE across all voxels', ...
        sprintf('TFCE FWE (Positive, %s)', tag), 'right', false);
    [hTfceNeg(side), results.pOmnibusNeg{side}] = ea_unified_perm_nulldist_plot( ...
        maxNegPerm, max(tfceEmpNeg, [], 'omitnan'), 'Max TFCE across all voxels (negated map)', ...
        sprintf('TFCE FWE (Negative, %s)', tag), 'right', false);

    fprintf('%s: %d/%d pos + %d/%d neg voxels survive TFCE-corrected p<=%.3g\n', ...
        tag, sum(survivePos), numel(survivePos), sum(surviveNeg), numel(surviveNeg), alpha);
end

ea_unified_perm_vals2nii(meta.space, tfceThreshMap, outdir, sprintf('%s_tfce_p%.3g', basefname, alpha));

for side = 1:nsides
    if nsides == 2
        sideTag = {'_r','_l'};
        sideTag = sideTag{side};
    else
        sideTag = '';
    end
    saveas(hTfcePos(side), fullfile(outdir, sprintf('%s_tfce_pos%s.png', basefname, sideTag)));
    saveas(hTfceNeg(side), fullfile(outdir, sprintf('%s_tfce_neg%s.png', basefname, sideTag)));
end

readmeBody = sprintf([ ...
    'TFCE (Threshold-Free Cluster Enhancement) analysis for %s.\n\n', ...
    'Uses PALM''s TFCE method (ea_unified_perm_palm_tfce.m, a port of PALM''s palm_tfce.m).\n', ...
    'If reporting this, cite:\n', ...
    '  Winkler AM, Ridgway GR, Douaud G, Nichols TE, Smith SM (2016). Faster permutation\n', ...
    '    inference in brain imaging. NeuroImage 141:502-516. (PALM)\n', ...
    '  Smith SM, Nichols TE (2009). Threshold-free cluster enhancement: addressing\n', ...
    '    problems of smoothing, threshold dependence and localisation in cluster\n', ...
    '    inference. NeuroImage 44:83-98. (the TFCE method itself)\n\n', ...
    'Settings used for this run:\n', ...
    '  E              = %g\n', ...
    '  H              = %g\n', ...
    '  conn           = %d\n', ...
    '  deltah         = %g (0 = auto, dh = max(map)/100)\n', ...
    '  NpermTFCE      = %d\n', ...
    '  alpha          = %.3g\n\n', ...
    'Files added by this step:\n\n', ...
    '  *_tfce_p%.3g_r/_l.nii\n', ...
    '    TFCE-enhanced map, thresholded at the per-voxel p-value computed by ranking\n', ...
    '    each voxel''s empirical TFCE value against the SAME null distribution (the\n', ...
    '    largest TFCE value anywhere in the map, one draw per permutation). THIS MAP\n', ...
    '    IS CORRECTED FOR MULTIPLE COMPARISONS (unlike *_voxelwise_uncorrected_* from\n', ...
    '    ea_unified_perm_sweetspot_stats -- it is a complementary alternative to that\n', ...
    '    file''s *_maxstat_* map, using cluster-aware enhancement instead of the raw\n', ...
    '    per-voxel statistic).\n\n', ...
    '  *_tfce_pos.png / *_tfce_neg.png\n', ...
    '    Null distribution of the largest TFCE value anywhere in the map, per\n', ...
    '    permutation, with this run''s own empirical whole-map max TFCE marked --\n', ...
    '    an omnibus "is there at least one significant TFCE-enhanced region anywhere"\n', ...
    '    readout, positive/negative directions plotted separately.\n'], ...
    basefname, E, H, conn, deltah, NpermTFCE, alpha, alpha);
ea_unified_perm_readme_append(outdir, 'TFCE analysis (ea_unified_perm_sweetspot_tfce)', readmeBody);

outfile = fullfile(outdir, 'tfce_results.mat');
results.savedfile = outfile;
save(outfile, '-struct', 'results', '-v7.3');
fprintf('Saved TFCE results to %s\n', outfile);
