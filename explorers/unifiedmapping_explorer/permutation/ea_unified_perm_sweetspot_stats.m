function results = ea_unified_perm_sweetspot_stats(permtestFile, alphaVoxelwise, alphaEisenstein, alphaMaxstat, useTailApprox, chunkSize)
% Reads a results.mat file (produced by ea_unified_perm_sweetspot) and
% computes three complementary permutation-based significance analyses,
% entirely standalone -- only the file path and explicit alpha levels are
% needed, no live ea_unifiedmapping object or GUI settings.
%
% Processes one side at a time, and within a side reads Rperm/pperm in
% voxel-column chunks (matfile row/column-chunked reads directly off disk,
% never materializing the full Nperm x Nvoxels arrays) -- this is the fix
% for a known memory problem in the ported-from version of this function,
% which loaded both sides' full Rperm/pperm at once and pushed MATLAB to
% ~159GB in testing at Nperm=5000 on a ~769K-voxel bilateral grid. Peak
% memory here is roughly Nperm x chunkSize regardless of Nperm/Nvoxels, a
% few hundred MB at the defaults below. The math/output is unchanged --
% only how it's computed changed. The trickiest part is the max-statistic
% test, which needs the max/min ACROSS ALL VOXELS per permutation: this is
% accumulated incrementally as chunks stream in (running max/min vector,
% O(Nperm) memory) rather than reduced in one shot.
%
% alphaVoxelwise  - per-voxel uncorrected permutation threshold. THIS TEST
%                    IS NOT CORRECTED FOR MULTIPLE COMPARISONS: each voxel
%                    is ranked only against its own, independent null
%                    distribution, with no accounting for the fact that
%                    ~769K such tests are being run simultaneously. At
%                    alpha=0.05 that scale, expect ~38K voxels to survive
%                    by chance alone even under a true global null (the
%                    same failure mode behind fMRI's well-known "dead
%                    salmon" cautionary result, Bennett et al. 2009, from
%                    skipping this exact correction). Treat the exported
%                    *_voxelwise_uncorrected_* map as EXPLORATORY/
%                    DESCRIPTIVE ONLY -- a picture of where the raw effect
%                    looks strongest, not a statistically defensible claim
%                    that any specific voxel is significant. The max-stat
%                    map below is the one that corrects for multiplicity.
%                    Default 0.05.
% alphaEisenstein - per-voxel significance cutoff used to decide which voxels
%                    are summed into the Eisenstein et al. 2014 (Annals of
%                    Neurology 76:279-295) omnibus Q statistic:
%                    Q = sum(-log10(p)) over voxels with p <= alphaEisenstein,
%                    computed for the empirical map and every permuted map,
%                    giving a permutation-based whole-map significance test
%                    (one test, not one per voxel -- properly controlled).
%                    Default 0.05.
% alphaMaxstat    - tail probability for the max-statistic FWER correction
%                    (e.g. 0.05 -> the 95th percentile of the per-permutation
%                    maximum, 5th percentile of the per-permutation minimum).
%                    This IS corrected for multiple comparisons: the null is
%                    built from the single most extreme voxel per
%                    permutation, the same principle FSL's randomise/PALM
%                    use. Default 0.05.
% useTailApprox   - if true, the Eisenstein and max-statistic p-values (each
%                    a single test against one shared null, unlike the
%                    per-voxel test) are refined with a Generalized Pareto
%                    Distribution tail fit when the rank-based p is already
%                    small (see ea_unified_perm_palm_gpd_pval, ported from FSL
%                    PALM's palm_pareto.m -- cite it if this changes a
%                    reported p-value). Default false.
% chunkSize       - voxels read per chunk. Default 20000 (a
%                    Nperm x chunkSize single-precision chunk is ~400MB at
%                    Nperm=5000); lower it further on tighter-memory
%                    machines, at the cost of more, smaller disk reads.
%
% Positive ("sweet spot") and negative ("sour spot") directions are tracked
% and thresholded separately throughout, rather than folded into |R|.
%
% Outputs, all written into the same folder permtestFile lives in (one
% *_r/_l nifti pair per map, one *_r/_l README entry per run):
%   *_voxelwise_uncorrected_p<alphaVoxelwise>_r/_l.nii - see alphaVoxelwise above.
%   *_maxstat_p<alphaMaxstat>_r/_l.nii                 - FWE-corrected map.
%   *_eisenstein.png                                    - omnibus null-distribution plot.
%   *_maxstat_pos.png / *_maxstat_neg.png               - max-stat null-distribution plots.
%   README.txt (appended)                               - settings used + what each file shows.
%
% Returns a struct with everything computed ({side}-indexed).

if ~exist('alphaVoxelwise', 'var') || isempty(alphaVoxelwise)
    alphaVoxelwise = 0.05;
end
if ~exist('alphaEisenstein', 'var') || isempty(alphaEisenstein)
    alphaEisenstein = 0.05;
end
if ~exist('alphaMaxstat', 'var') || isempty(alphaMaxstat)
    alphaMaxstat = 0.05;
end
if ~exist('useTailApprox', 'var') || isempty(useTailApprox)
    useTailApprox = false;
end
if ~exist('chunkSize', 'var') || isempty(chunkSize)
    chunkSize = 20000;
end

meta = load(permtestFile, 'space', 'Remp', 'pemp', 'Nperm', 'nsides');
mFile = matfile(permtestFile);
outdir = [fileparts(permtestFile), filesep];
[~, basefname] = fileparts(fileparts(permtestFile)); % permtestFile lives in its own run folder -- use that folder's name to tag exports

nsides = meta.nsides;
Nperm = meta.Nperm;

results = struct();
results.sourcefile = permtestFile;
results.alphaVoxelwise = alphaVoxelwise;
results.alphaEisenstein = alphaEisenstein;
results.alphaMaxstat = alphaMaxstat;
results.useTailApprox = useTailApprox;

voxThreshMap = cell(1, nsides);
maxstatMap = cell(1, nsides);
hEisenstein = gobjects(1, nsides);
hMaxstatPos = gobjects(1, nsides);
hMaxstatNeg = gobjects(1, nsides);

for side = 1:nsides
    Remp = meta.Remp{side};   % Nvoxels x 1
    pemp = meta.pemp{side};   % Nvoxels x 1
    Nvox = numel(Remp);
    RempRow = Remp'; % 1 x Nvoxels, for broadcasting against each chunk's rows

    RvarName = sprintf('Rperm_side%d', side);
    pvarName = sprintf('pperm_side%d', side);

    countPos = zeros(1, Nvox);
    countNeg = zeros(1, Nvox);
    Qperm = zeros(Nperm, 1);          % running sum across voxel chunks
    maxRperm = -inf(Nperm, 1);        % running max across voxel chunks
    minRperm = inf(Nperm, 1);         % running min across voxel chunks

    tag = sprintf('side%d', side);
    fprintf('%s: processing %d voxels in chunks of %d...\n', tag, Nvox, chunkSize);

    for vStart = 1:chunkSize:Nvox
        vEnd = min(vStart+chunkSize-1, Nvox);
        vIdx = vStart:vEnd;

        Rchunk = double(mFile.(RvarName)(:, vIdx)); % Nperm x numel(vIdx)
        pchunk = double(mFile.(pvarName)(:, vIdx));

        %% 1. Voxelwise uncorrected permutation threshold (per-voxel own null, pos/neg separate)
        countPos(vIdx) = countPos(vIdx) + sum(Rchunk >= RempRow(vIdx), 1);
        countNeg(vIdx) = countNeg(vIdx) + sum(Rchunk <= RempRow(vIdx), 1);

        %% 2. Eisenstein 2014 omnibus Q: running sum of -log10(p) over per-voxel-significant voxels
        sigMaskPerm = pchunk <= alphaEisenstein;
        logpPerm = -log10(max(pchunk, eps));
        logpPerm(~sigMaskPerm) = 0;
        Qperm = Qperm + ea_nansum(logpPerm, 2); % reduce over THIS CHUNK's voxels (dim 2), accumulate across chunks

        %% 3. Max-statistic FWER: running max/min across voxel chunks
        % ea_nanmax/ea_nanmin (ext_libs/nan) do NOT use MATLAB's own
        % max(A,[],dim) convention -- their 2-arg form (a,dim) reduces along
        % a dimension (used here to collapse this chunk's voxels), while
        % their 3-arg form with dim=[] is an ELEMENTWISE max/min of two
        % same-sized arrays (used here to combine the running accumulator
        % with this chunk's reduced result). Mixing these up has caused a
        % real crash elsewhere in this codebase -- see the two distinct
        % calls below, not interchangeable.
        chunkMax = ea_nanmax(Rchunk, 2); % Nperm x 1, reduce along dim 2 (voxels in this chunk)
        chunkMin = ea_nanmin(Rchunk, 2);
        maxRperm = ea_nanmax(maxRperm, [], chunkMax); % Nperm x 1, elementwise combine with running max
        minRperm = ea_nanmin(minRperm, [], chunkMin); % Nperm x 1, elementwise combine with running min
    end

    % +1 in numerator and denominator: the empirical map is itself one
    % valid draw under the null, so it belongs in its own reference set.
    % Without it, a voxel with zero exceedances would report p = 0, which
    % is not a valid claim from a finite number of permutations (Phipson &
    % Smyth 2010).
    pVoxPos = (countPos + 1) / (Nperm + 1);
    pVoxNeg = (countNeg + 1) / (Nperm + 1);
    pVoxPos(isnan(RempRow)) = NaN;
    pVoxNeg(isnan(RempRow)) = NaN;

    voxSurvivePos = pVoxPos <= alphaVoxelwise;
    voxSurviveNeg = pVoxNeg <= alphaVoxelwise;

    thisVoxMap = nan(size(Remp));
    thisVoxMap(voxSurvivePos') = Remp(voxSurvivePos');
    thisVoxMap(voxSurviveNeg') = Remp(voxSurviveNeg');
    voxThreshMap{side} = thisVoxMap;

    results.voxelwise.pPos{side} = pVoxPos;
    results.voxelwise.pNeg{side} = pVoxNeg;
    results.voxelwise.thresholdedMap{side} = thisVoxMap;

    %% Eisenstein: empirical Q, then rank Qperm (already fully reduced -- small, Nperm x 1) against it
    sigMaskEmp = pemp' <= alphaEisenstein;
    logpEmp = -log10(max(pemp', eps));
    logpEmp(~sigMaskEmp) = 0;
    Qemp = ea_nansum(logpEmp);

    results.eisenstein.Qemp{side} = Qemp;
    results.eisenstein.Qperm{side} = Qperm;

    [hEisenstein(side), results.eisenstein.p{side}, results.eisenstein.rank{side}] = ea_unified_perm_nulldist_plot( ...
        Qperm, Qemp, 'Q = \Sigma -log_{10}(p)', ...
        sprintf('Eisenstein 2014 Omnibus Test (%s, \\alpha=%.3g)', tag, alphaEisenstein), 'right', useTailApprox);

    %% Max-statistic: empirical max/min, then rank the (already fully accumulated) null against it
    maxRemp = ea_nanmax(Remp);
    minRemp = ea_nanmin(Remp);

    % Exact-rank threshold, consistent with the (count+1)/(Nperm+1)
    % p-value formula used above and in ea_unified_perm_nulldist_plot --
    % NOT prctile(), whose interpolated quantile can land on a value that
    % never actually occurred in the null, so the exported thresholded map
    % and the printed p-value for the same test could silently disagree at
    % the edge.
    threshPos = ea_unified_perm_exact_maxstat_threshold(maxRperm, alphaMaxstat, 'right');
    threshNeg = ea_unified_perm_exact_maxstat_threshold(minRperm, alphaMaxstat, 'left');

    thisMaxstatMap = nan(size(Remp));
    thisMaxstatMap(Remp >= threshPos) = Remp(Remp >= threshPos);
    thisMaxstatMap(Remp <= threshNeg) = Remp(Remp <= threshNeg);
    maxstatMap{side} = thisMaxstatMap;

    results.maxstat.threshPos{side} = threshPos;
    results.maxstat.threshNeg{side} = threshNeg;
    results.maxstat.thresholdedMap{side} = thisMaxstatMap;

    [hMaxstatPos(side), results.maxstat.pPos{side}, results.maxstat.rankPos{side}] = ea_unified_perm_nulldist_plot( ...
        maxRperm, maxRemp, 'Max R across all voxels', ...
        sprintf('Max-Statistic FWER (Positive, %s)', tag), 'right', useTailApprox);

    [hMaxstatNeg(side), results.maxstat.pNeg{side}, results.maxstat.rankNeg{side}] = ea_unified_perm_nulldist_plot( ...
        minRperm, minRemp, 'Min R across all voxels', ...
        sprintf('Max-Statistic FWER (Negative, %s)', tag), 'left', useTailApprox);

    fprintf('%s: voxelwise (UNCORRECTED) %d/%d pos + %d/%d neg voxels survive p<=%.3g | Eisenstein Q p=%.3g | max-stat (FWE-corrected) pos p=%.3g, neg p=%.3g\n', ...
        tag, sum(voxSurvivePos), numel(voxSurvivePos), sum(voxSurviveNeg), numel(voxSurviveNeg), alphaVoxelwise, ...
        results.eisenstein.p{side}, results.maxstat.pPos{side}, results.maxstat.pNeg{side});
end

ea_unified_perm_vals2nii(meta.space, voxThreshMap, outdir, sprintf('%s_voxelwise_uncorrected_p%.3g', basefname, alphaVoxelwise));
ea_unified_perm_vals2nii(meta.space, maxstatMap, outdir, sprintf('%s_maxstat_p%.3g', basefname, alphaMaxstat));

% Save the null-distribution figures alongside everything else in this
% run's folder (they're otherwise only ever displayed, never persisted) --
% one file per side if there are two.
for side = 1:nsides
    if nsides == 2
        sideTag = {'_r','_l'};
        sideTag = sideTag{side};
    else
        sideTag = '';
    end
    saveas(hEisenstein(side), fullfile(outdir, sprintf('%s_eisenstein%s.png', basefname, sideTag)));
    saveas(hMaxstatPos(side), fullfile(outdir, sprintf('%s_maxstat_pos%s.png', basefname, sideTag)));
    saveas(hMaxstatNeg(side), fullfile(outdir, sprintf('%s_maxstat_neg%s.png', basefname, sideTag)));
end

readmeBody = sprintf([ ...
    'Post-hoc statistical report for %s (source: %s).\n\n', ...
    'Settings used for this run:\n', ...
    '  alphaVoxelwise  = %.3g\n', ...
    '  alphaEisenstein = %.3g\n', ...
    '  alphaMaxstat    = %.3g\n', ...
    '  useTailApprox   = %d\n\n', ...
    'Files added by this step:\n\n', ...
    '  *_voxelwise_uncorrected_p%.3g_r/_l.nii\n', ...
    '    Per-voxel permutation test: each voxel''s real R value is compared ONLY\n', ...
    '    against its own null distribution (that voxel''s Nperm permuted R values),\n', ...
    '    separately for positive (''sweet spot'': more field -> better outcome) and\n', ...
    '    negative (''sour spot'': more field -> worse outcome) directions.\n', ...
    '    SURVIVAL CONDITION for a voxel to appear non-NaN in this map:\n', ...
    '      pPos = (count of permuted R >= real R, at that voxel, +1) / (Nperm+1) <= alphaVoxelwise\n', ...
    '      -- OR --\n', ...
    '      pNeg = (count of permuted R <= real R, at that voxel, +1) / (Nperm+1) <= alphaVoxelwise\n', ...
    '    THIS TEST IS NOT CORRECTED FOR MULTIPLE COMPARISONS. Every voxel is its own\n', ...
    '    independent test, with no accounting for the fact that hundreds of thousands\n', ...
    '    of such tests are run at once here. At alpha=%.3g across a ~769K-voxel grid,\n', ...
    '    expect on the order of tens of thousands of voxels to survive by chance alone\n', ...
    '    even if there is truly no effect anywhere. This map is EXPLORATORY/\n', ...
    '    DESCRIPTIVE ONLY -- a picture of where the raw effect looks strongest -- and\n', ...
    '    should NOT be reported as a standalone claim that any specific voxel is\n', ...
    '    statistically significant. Use the maxstat map below for that claim.\n\n', ...
    '  *_maxstat_p%.3g_r/_l.nii\n', ...
    '    Max-statistic family-wise error (FWE) corrected map: the null distribution\n', ...
    '    here is the SINGLE most extreme voxel value, per permutation, across the\n', ...
    '    WHOLE map (not per-voxel) -- the same principle used by FSL''s randomise/\n', ...
    '    PALM. This DOES correct for multiple comparisons. A voxel survives if its\n', ...
    '    real R exceeds the alphaMaxstat-level threshold of that whole-map null\n', ...
    '    (separately for the positive/max and negative/min tails). This is the map\n', ...
    '    that supports a statistically defensible significance claim.\n\n', ...
    '  *_eisenstein.png\n', ...
    '    Null distribution for the Eisenstein et al. 2014 omnibus test: a SINGLE\n', ...
    '    whole-map statistic, Q = sum(-log10(p)) over voxels with p<=alphaEisenstein,\n', ...
    '    computed once for the real map and once per permutation. Answers "does this\n', ...
    '    map contain more low-p voxels overall than chance would produce" -- a\n', ...
    '    properly-corrected, one-test-total alternative to reading the uncorrected\n', ...
    '    voxelwise map. Does NOT tell you which specific voxels are responsible.\n\n', ...
    '  *_maxstat_pos.png / *_maxstat_neg.png\n', ...
    '    Null distributions underlying the maxstat nifti above, positive/negative\n', ...
    '    tails plotted separately, empirical value marked.\n'], ...
    basefname, permtestFile, alphaVoxelwise, alphaEisenstein, alphaMaxstat, useTailApprox, ...
    alphaVoxelwise, alphaVoxelwise, alphaMaxstat);
ea_unified_perm_readme_append(outdir, 'Post-hoc statistics (ea_unified_perm_sweetspot_stats)', readmeBody);
end

function thresh = ea_unified_perm_exact_maxstat_threshold(nullvals, alpha, tail)
% Exact-rank threshold matching the (count+1)/(Nperm+1) p-value formula used
% throughout this file and in ea_unified_perm_nulldist_plot -- the value at
% which a tied empirical statistic would first report p <= alpha, using
% only order statistics actually present in nullvals (unlike prctile, which
% interpolates between them).
%
% tail: 'right' (large nullvals are extreme, e.g. positive max-stat)
%       'left'  (small nullvals are extreme, e.g. negative max-stat)
%
% Returns NaN if no achievable rank gets p <= alpha with this many
% permutations (i.e. even the single most extreme null value isn't
% enough) -- comparisons against NaN are always false, so downstream
% thresholding correctly marks nothing as surviving.

Nperm = numel(nullvals);
switch tail
    case 'right'
        sorted = sort(nullvals, 'descend');
    case 'left'
        sorted = sort(nullvals, 'ascend');
end

pAtRank = ((1:Nperm) + 1) / (Nperm + 1); % p-value if the empirical stat ties the k-th most extreme null value
survivingIdx = find(pAtRank <= alpha, 1, 'last');
if isempty(survivingIdx)
    thresh = NaN;
else
    thresh = sorted(survivingIdx);
end
end
