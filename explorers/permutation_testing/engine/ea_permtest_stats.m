function stats = ea_permtest_stats(mapset, opts)
% Permutation statistics on a "map set" (see ea_permtest_read_fiberrun.m):
% shared by all modalities -- the units of a map can be fibers or voxels,
% the calculations are identical. Large shuffled matrices are read from
% disk in chunks of units, so memory use does not grow with map size.
%
% 1. Eisenstein omnibus test (Eisenstein et al. 2014, Ann Neurol 76:279-295):
%    one test per map/analysis: "does this map contain more association
%    with outcome than chance?". Per unit the parametric p of the
%    statistic is used:
%      Q    = sum of -log10(p)   over units with p <= alphaEisenstein
%      Qpos = sum of -log10(p/2) over units with R > 0 and p/2 <= alphaEisenstein
%             (one-sided, improvement direction)
%      Qneg = same for R < 0
%    computed for the real map and for every shuffled map; omnibus
%    p = (#shuffles with Q >= real Q + 1) / (Nperm + 1).
% 2. Max statistic (Nichols & Holmes 2002), family-wise error (FWE)
%    corrected per unit: the null distribution is the highest statistic
%    of each shuffled map over all units of the family; a unit survives
%    if its real statistic beats (1 - alpha) of these maxima,
%    p_FWE = (#maxima >= real + 1) / (Nperm + 1).
%    - single-step: one null per family;
%    - step-down (Westfall & Young 1993): units already found significant
%      are removed and the maxima recomputed over the remaining units,
%      repeated until nothing new survives. Same FWE control, never fewer
%      significant units. Implemented iteratively (one pass over the data
%      per iteration) so it also works for voxel maps that do not fit in
%      memory. The set of significant units is exact; their step-down
%      p-values are exact (successive maxima) when their shuffled values
%      fit in opts.maxExactStepdownMB, otherwise conservative; for
%      non-significant units the reported step-down p is a conservative
%      upper bound (always > alpha).
%    Families:
%      'perside'     each map separately (e.g. one subscore, one side)
%      'bothsides'   all maps of one analysis (e.g. one subscore, L + R)
%      'allanalyses' all maps (e.g. all subscores, L + R)
%
% Rank-based p-values never reach 0 (Phipson & Smyth 2010); a shuffle
% without any defined statistic never counts as exceeding but stays in
% the denominator (same conventions as the sweetspot permutation code).
%
% mapset - from a reader (ea_permtest_read_*.m).
% opts   - struct, all optional:
%   .alpha            FWE level for the max statistic. Default 0.05.
%   .alphaEisenstein  per-unit p cutoff entering Q. Default 0.05.
%   .direction        'positive' (default; R > 0 only, at alpha),
%                     'negative', 'both-each' (each direction at alpha, as
%                     in the sweetspot statistics), 'both-split' (each at
%                     alpha/2, 5% overall).
%   .families         cellstr, default {'perside','bothsides','allanalyses'}.
%   .omnibus          cellstr of 'Q','Qpos','Qneg', default {'Q','Qpos'}.
%   .omnibusScopes    cellstr of 'perside','bothsides', default both.
%   .stepdown         default true.
%   .useTailApprox    GPD tail fit for omnibus and single-step p-values
%                     (ea_permtest_gpd_pval). Default false.
%   .chunkSize        units read per chunk. Default 5000.
%   .maxExactStepdownMB  memory allowed for exact step-down p-values of
%                     significant units. Default 2000.
%   .plots            write PNG figures. Default true.
%   .outdir           output folder (created). Default: a new subfolder of
%                     mapset.statsRoot named after the settings and time.
%
% stats - struct with .settings, .maps (Remp per unit), .omnibus(k) and
%         .maxstat(k) (see the fields set below); saved as stats.mat in
%         outdir, summarised in the run's README.md.

if ~exist('opts', 'var') || isempty(opts)
    opts = struct;
end
opts = local_defaults(opts, struct('alpha', 0.05, 'alphaEisenstein', 0.05, 'direction', 'positive', ...
    'families', {{'perside', 'bothsides', 'allanalyses'}}, 'omnibus', {{'Q', 'Qpos'}}, ...
    'omnibusScopes', {{'perside', 'bothsides'}}, 'stepdown', true, 'useTailApprox', false, ...
    'chunkSize', 5000, 'maxExactStepdownMB', 2000, 'plots', true, 'outdir', ''));

switch opts.direction
    case 'positive'
        tails = {'pos'}; alphaTail = opts.alpha;
    case 'negative'
        tails = {'neg'}; alphaTail = opts.alpha;
    case 'both-each'
        tails = {'pos', 'neg'}; alphaTail = opts.alpha;
    case 'both-split'
        tails = {'pos', 'neg'}; alphaTail = opts.alpha / 2;
    otherwise
        ea_error('opts.direction must be ''positive'', ''negative'', ''both-each'' or ''both-split''.');
end

maps = mapset.maps;
Nmaps = numel(maps);
if Nmaps == 0
    ea_error('The map set contains no units to test.');
end

if isempty(opts.outdir)
    opts.outdir = fullfile(mapset.statsRoot, sprintf('%s_alpha%g_%s', opts.direction, opts.alpha, ...
        datestr(now, 'yyyymmdd_HHMMSS')));
end
[ok, msg] = mkdir(opts.outdir);
if ~ok
    ea_error(['Cannot create the statistics folder:\n', local_esc(opts.outdir), '\n', local_esc(msg)]);
end
fprintf('Permutation statistics -> %s\n', opts.outdir);

%% Pass 1: per map, per-shuffle maxima (both tails) and omnibus Q of every shuffle
fprintf('Pass 1: maxima and omnibus statistics of every shuffle...\n');
mapNull = struct('pos', cell(1, Nmaps), 'neg', []);
mapQ = struct('perm', cell(1, Nmaps), 'emp', []);
for i = 1:Nmaps
    N = maps(i).Nperm;
    nullPos = -inf(N, 1);
    nullNeg = -inf(N, 1);
    Qperm = struct();
    for t = opts.omnibus
        Qperm.(t{1}) = zeros(N, 1);
    end
    m = matfile(maps(i).file);
    for c0 = 1:opts.chunkSize:maps(i).Nunits
        cols = c0:min(c0 + opts.chunkSize - 1, maps(i).Nunits);
        R = m.(maps(i).Rvar)(:, cols);
        nullPos = max(nullPos, max(R, [], 2));
        nullNeg = max(nullNeg, max(-R, [], 2));
        if ~isempty(opts.omnibus)
            P = double(m.(maps(i).pvar)(:, cols));
            for t = opts.omnibus
                Qperm.(t{1}) = Qperm.(t{1}) + local_Q(t{1}, double(R), P, opts.alphaEisenstein);
            end
        end
    end
    mapNull(i).pos = nullPos;
    mapNull(i).neg = nullNeg;
    mapQ(i).perm = Qperm;
    Qemp = struct();
    for t = opts.omnibus
        Qemp.(t{1}) = local_Q(t{1}, maps(i).Remp(:)', maps(i).pemp(:)', opts.alphaEisenstein);
    end
    mapQ(i).emp = Qemp;
    fprintf('  %s: %d units\n', maps(i).label, maps(i).Nunits);
end

%% Omnibus tests
omnibus = struct('type', {}, 'scope', {}, 'analysis', {}, 'members', {}, 'label', {}, ...
    'emp', {}, 'null', {}, 'p', {}, 'tailFitUsed', {});
analysisIdx = [maps.analysis];
Nanalyses = numel(mapset.analyses);
for a = 1:Nanalyses
    members = find(analysisIdx == a);
    if isempty(members)
        continue
    end
    for t = opts.omnibus
        if ismember('bothsides', opts.omnibusScopes) && local_sameN(maps(members))
            nullQ = zeros(maps(members(1)).Nperm, 1);
            empQ = 0;
            for i = members
                nullQ = nullQ + mapQ(i).perm.(t{1});
                empQ = empQ + mapQ(i).emp.(t{1});
            end
            omnibus(end+1) = local_omnibus(t{1}, 'bothsides', a, members, ...
                sprintf('%s, %s', mapset.analyses(a).label, local_scopeLabel(maps(members))), empQ, nullQ, opts); %#ok<AGROW>
        end
        if ismember('perside', opts.omnibusScopes)
            for i = members
                omnibus(end+1) = local_omnibus(t{1}, 'perside', a, i, maps(i).label, ...
                    mapQ(i).emp.(t{1}), mapQ(i).perm.(t{1}), opts); %#ok<AGROW>
            end
        end
    end
end

%% Max statistic: families
fams = struct('family', {}, 'members', {}, 'label', {}, 'analysis', {});
if ismember('perside', opts.families)
    for i = 1:Nmaps
        fams(end+1) = struct('family', 'perside', 'members', i, 'label', maps(i).label, 'analysis', maps(i).analysis); %#ok<AGROW>
    end
end
if ismember('bothsides', opts.families)
    for a = 1:Nanalyses
        members = find(analysisIdx == a);
        if ~isempty(members) && local_sameN(maps(members))
            fams(end+1) = struct('family', 'bothsides', 'members', members, ...
                'label', sprintf('%s, %s', mapset.analyses(a).label, local_scopeLabel(maps(members))), 'analysis', a); %#ok<AGROW>
        end
    end
end
if ismember('allanalyses', opts.families)
    if local_sameN(maps)
        fams(end+1) = struct('family', 'allanalyses', 'members', 1:Nmaps, ...
            'label', sprintf('all %d analyses, all sides', Nanalyses), 'analysis', 0);
    else
        ea_cprintf('CmdWinWarnings', ['Skipping the ''allanalyses'' family: the analyses have different ', ...
            'numbers of shuffles (small groups enumerated exactly).\n']);
    end
end

maxstat = struct('family', {}, 'tail', {}, 'alpha', {}, 'members', {}, 'label', {}, 'analysis', {}, ...
    'Nperm', {}, 'null', {}, 'threshold', {}, 'observedMax', {}, 'familyP', {}, 'familyPTailFit', {}, ...
    'pSingle', {}, 'sigSingle', {}, 'NsigSingle', {}, 'pStepdown', {}, 'sigStepdown', {}, ...
    'NsigStepdown', {}, 'stepdownIterations', {}, 'stepdownExact', {});
for t = tails
    sgn = 1 - 2 * strcmp(t{1}, 'neg'); % +1 for pos, -1 for neg
    for k = 1:numel(fams)
        fam = fams(k);
        members = fam.members;
        N = maps(members(1)).Nperm;
        nullv = -inf(N, 1);
        for i = members
            nullv = max(nullv, mapNull(i).(t{1}));
        end
        x = cell(1, numel(members));
        for j = 1:numel(members)
            x{j} = double(single(sgn * maps(members(j)).Remp(:))); % compared at the stored (single) precision
        end
        [pS, sigS] = local_pvals(nullv, x, alphaTail, opts);
        thr = local_threshold(nullv, alphaTail);
        allx = vertcat(x{:});
        [obsMax, ~] = max(allx);
        [famP, famTail] = local_rankp(nullv, obsMax, opts.useTailApprox);

        e = struct('family', fam.family, 'tail', t{1}, 'alpha', alphaTail, 'members', members, ...
            'label', fam.label, 'analysis', fam.analysis, 'Nperm', N, 'null', sgn * nullv, ...
            'threshold', sgn * thr, 'observedMax', sgn * obsMax, 'familyP', famP, 'familyPTailFit', famTail, ...
            'pSingle', {pS}, 'sigSingle', {sigS}, 'NsigSingle', sum(cellfun(@nnz, sigS)), ...
            'pStepdown', {pS}, 'sigStepdown', {sigS}, 'NsigStepdown', sum(cellfun(@nnz, sigS)), ...
            'stepdownIterations', 0, 'stepdownExact', true);
        if opts.stepdown
            [e.pStepdown, e.sigStepdown, e.stepdownIterations, e.stepdownExact] = ...
                local_stepdown(maps, members, sgn, x, sigS, pS, alphaTail, opts);
            e.NsigStepdown = sum(cellfun(@nnz, e.sigStepdown));
        end
        maxstat(end+1) = e; %#ok<AGROW>
        fprintf('  max-stat %-9s %s %-40s: %d significant (single-step), %d (step-down), threshold R = %.3f\n', ...
            fam.family, t{1}, fam.label, e.NsigSingle, e.NsigStepdown, e.threshold);
    end
end

%% Save, plot, document
stats = struct;
stats.settings = opts;
stats.settings.tails = tails;
stats.settings.alphaPerTail = alphaTail;
stats.modality = mapset.modality;
stats.source = mapset.source;
stats.analyses = mapset.analyses;
stats.maps = rmfield(maps, {'pemp'});
stats.omnibus = omnibus;
stats.maxstat = maxstat;
stats.created = datestr(now, 'yyyy-mm-dd HH:MM:SS');
save(fullfile(opts.outdir, 'stats.mat'), 'stats', '-v7.3');

if opts.plots
    local_plots(stats, mapset, opts);
end
local_readme(stats, mapset, opts);
fprintf('Statistics saved in %s\n', opts.outdir);


%% ------------------------------------------------------------------------
function Q = local_Q(type, R, P, alphaE)
% Omnibus Q per row (shuffle) for a block of units. R, P: rows x units.
switch type
    case 'Q'
        Pt = P;
        keep = Pt <= alphaE;
    case 'Qpos'
        Pt = P / 2;
        keep = R > 0 & Pt <= alphaE;
    case 'Qneg'
        Pt = P / 2;
        keep = R < 0 & Pt <= alphaE;
    otherwise
        ea_error(['Unknown omnibus type ', type]);
end
L = -log10(max(Pt, eps));
L(~keep) = 0;
Q = sum(L, 2);


function o = local_omnibus(type, scope, a, members, label, empQ, nullQ, opts)
[p, tailFit] = local_rankp(nullQ, empQ, opts.useTailApprox);
o = struct('type', type, 'scope', scope, 'analysis', a, 'members', members, 'label', label, ...
    'emp', empQ, 'null', nullQ, 'p', p, 'tailFitUsed', tailFit);


function [p, tailFit] = local_rankp(nullv, empv, useTailApprox)
p = (sum(nullv >= empv) + 1) / (numel(nullv) + 1);
tailFit = false;
if useTailApprox && isfinite(empv)
    [p, tailFit] = ea_permtest_gpd_pval(nullv, empv, 'right');
end


function [pc, sigc] = local_pvals(nullv, x, alpha, opts)
% Single-step FWE p per unit against one null (rank-based, optional GPD
% refinement of small p-values).
N = numel(nullv);
pc = cell(1, numel(x));
sigc = cell(1, numel(x));
for j = 1:numel(x)
    xj = x{j};
    p = nan(numel(xj), 1);
    for c0 = 1:20000:numel(xj)
        cols = c0:min(c0 + 19999, numel(xj));
        cnt = sum(nullv >= xj(cols)', 1)';
        p(cols) = (cnt + 1) / (N + 1);
    end
    p(isnan(xj)) = NaN;
    if opts.useTailApprox
        for u = find(p < 0.10)'
            p(u) = ea_permtest_gpd_pval(nullv, xj(u), 'right');
        end
    end
    pc{j} = p;
    sigc{j} = p <= alpha;
end


function thr = local_threshold(nullv, alpha)
% Smallest real value that is still NOT significant: a unit is significant
% iff its value is above thr (exact rank, consistent with (count+1)/(N+1)).
N = numel(nullv);
kAllowed = floor(alpha * (N + 1) - 1 + 1e-9); % maxima allowed to reach the value
s = sort(nullv, 'descend');
if kAllowed < 0
    thr = inf;
elseif kAllowed + 1 > N
    thr = -inf;
else
    thr = s(kAllowed + 1);
end


function [pD, sigD, nIter, exact] = local_stepdown(maps, members, sgn, x, sigS, pS, alpha, opts)
% Iterative Westfall-Young step-down (see header).
nm = numel(members);
S = sigS;
pD = pS;
sigD = sigS;
nIter = 1;
exact = true;
if ~any(cellfun(@any, S))
    return % nothing significant: step-down equals single-step
end
while true
    nIter = nIter + 1;
    [nullK, ~] = local_pass(maps, members, sgn, S, false, opts);
    newAny = false;
    for j = 1:nm
        rest = ~S{j} & ~isnan(x{j});
        cnt = zeros(nnz(rest), 1);
        xr = x{j}(rest);
        for c0 = 1:20000:numel(xr)
            cols = c0:min(c0 + 19999, numel(xr));
            cnt(cols) = sum(nullK >= xr(cols)', 1)';
        end
        p = (cnt + 1) / (numel(nullK) + 1);
        pD{j}(rest) = p;
        newSig = false(size(S{j}));
        newSig(rest) = p <= alpha;
        if any(newSig)
            newAny = true;
            S{j} = S{j} | newSig;
        end
    end
    if ~newAny
        break
    end
end

% Exact successive-maxima p-values for the significant units, if their
% shuffled values fit in memory
nS = sum(cellfun(@nnz, S));
Nperm = maps(members(1)).Nperm;
if nS * Nperm * 4 / 2^20 <= opts.maxExactStepdownMB
    [nullRest, C] = local_pass(maps, members, sgn, S, true, opts);
    xs = [];
    for j = 1:nm
        xs = [xs; x{j}(S{j})]; %#ok<AGROW>
    end
    Cm = [C{:}];
    Cm(isnan(Cm)) = -inf; % an undefined shuffled value never exceeds
    [xsSorted, ord] = sort(xs, 'descend');
    Cm = Cm(:, ord);
    succMax = max(cummax(Cm, 2, 'reverse'), nullRest); % max over the unit and all lower-ranked units
    pExact = (sum(succMax >= xsSorted', 1)' + 1) / (Nperm + 1);
    pExact = cummax(pExact); % monotone in rank
    pBack = zeros(size(pExact));
    pBack(ord) = pExact;
    pos = 0;
    for j = 1:nm
        n = nnz(S{j});
        pD{j}(S{j}) = pBack(pos + (1:n));
        pos = pos + n;
    end
else
    exact = false; % iteration p-values: valid but conservative
end
% Monotonicity: no non-significant unit gets a smaller p than a significant one
pmaxS = max(cell2mat(cellfun(@(p, s) p(s), pD, S, 'Uni', 0)'));
for j = 1:nm
    rest = ~S{j} & ~isnan(x{j});
    pD{j}(rest) = max(pD{j}(rest), pmaxS);
end
sigD = cellfun(@(p) p <= alpha, pD, 'Uni', 0);


function [nullv, C] = local_pass(maps, members, sgn, S, collect, opts)
% One chunked pass over the family: per-shuffle maximum over units NOT in
% S; optionally collect the shuffled values of the units in S.
N = maps(members(1)).Nperm;
nullv = -inf(N, 1);
C = cell(1, numel(members));
for j = 1:numel(members)
    i = members(j);
    m = matfile(maps(i).file);
    Cj = {};
    for c0 = 1:opts.chunkSize:maps(i).Nunits
        cols = c0:min(c0 + opts.chunkSize - 1, maps(i).Nunits);
        R = sgn * m.(maps(i).Rvar)(:, cols);
        ex = S{j}(cols);
        if any(ex)
            if collect
                Cj{end+1} = double(R(:, ex)); %#ok<AGROW>
            end
            R(:, ex) = -inf;
        end
        nullv = max(nullv, max(R, [], 2));
    end
    if collect
        if isempty(Cj)
            C{j} = zeros(N, 0);
        else
            C{j} = [Cj{:}];
        end
    end
end
nullv = double(nullv);


function ok = local_sameN(mapsSub)
ok = numel(unique([mapsSub.Nperm])) == 1;


function s = local_scopeLabel(mapsSub)
if numel(mapsSub) > 1
    s = 'both sides';
else
    s = mapsSub(1).sideLabel;
end


function local_plots(stats, mapset, opts)
pal = ea_permtest_palette(); %#ok<NASGU>
for a = 1:numel(mapset.analyses)
    fname = regexprep(mapset.analyses(a).label, '[^\w-]', '_');
    % Omnibus: rows = Q types, columns = scopes
    om = stats.omnibus([stats.omnibus.analysis] == a);
    if ~isempty(om)
        types = unique({om.type}, 'stable');
        scopes = unique({om.label}, 'stable');
        fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 420 * numel(scopes), 320 * numel(types)]);
        tl = tiledlayout(fig, numel(types), numel(scopes), 'TileSpacing', 'compact', 'Padding', 'compact');
        for r = 1:numel(types)
            for c = 1:numel(scopes)
                o = om(strcmp({om.type}, types{r}) & strcmp({om.label}, scopes{c}));
                ax = nexttile(tl);
                if isempty(o)
                    axis(ax, 'off');
                    continue
                end
                ea_permtest_nulldist_panel(ax, o.null, o.emp, sprintf('%s = \\Sigma -log_{10}(p)', o.type), ...
                    sprintf('%s: %s', o.type, o.label), opts.useTailApprox);
            end
        end
        title(tl, sprintf('Omnibus test (Eisenstein 2014), \\alpha_{unit} = %g', opts.alphaEisenstein));
        exportgraphics(fig, fullfile(opts.outdir, sprintf('%s_omnibus.png', fname)), 'Resolution', 300);
        close(fig);
    end
    % Max statistic: rows = tails, columns = families of this analysis
    ms = stats.maxstat([stats.maxstat.analysis] == a);
    if ~isempty(ms)
        local_maxstat_figure(ms, fullfile(opts.outdir, sprintf('%s_maxstat.png', fname)), opts);
    end
end
ms = stats.maxstat(strcmp({stats.maxstat.family}, 'allanalyses'));
if ~isempty(ms)
    local_maxstat_figure(ms, fullfile(opts.outdir, 'allanalyses_maxstat.png'), opts);
end


function local_maxstat_figure(ms, file, opts)
tails = unique({ms.tail}, 'stable');
labels = unique({ms.label}, 'stable');
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 420 * numel(labels), 320 * numel(tails)]);
tl = tiledlayout(fig, numel(tails), numel(labels), 'TileSpacing', 'compact', 'Padding', 'compact');
for r = 1:numel(tails)
    for c = 1:numel(labels)
        e = ms(strcmp({ms.tail}, tails{r}) & strcmp({ms.label}, labels{c}));
        ax = nexttile(tl);
        if isempty(e)
            axis(ax, 'off');
            continue
        end
        sgn = 1 - 2 * strcmp(e.tail, 'neg');
        if sgn > 0
            xl = 'highest R per shuffle';
        else
            xl = '-(lowest R) per shuffle';
        end
        ea_permtest_nulldist_panel(ax, sgn * e.null, sgn * e.observedMax, xl, ...
            sprintf('%s: %d significant, %d with step-down', e.label, e.NsigSingle, e.NsigStepdown), ...
            opts.useTailApprox, sgn * e.threshold, sprintf('FWE threshold, \\alpha = %g', e.alpha));
    end
end
title(tl, 'Max statistic (FWE-corrected)');
exportgraphics(fig, file, 'Resolution', 300);
close(fig);


function local_readme(stats, mapset, opts)
L = {};
L{end+1} = sprintf('- Folder: `%s`', opts.outdir);
L{end+1} = sprintf('- Direction: %s, FWE alpha per direction: %g, omnibus unit cutoff p <= %g, step-down: %d, GPD tail fit: %d', ...
    opts.direction, stats.settings.alphaPerTail, opts.alphaEisenstein, opts.stepdown, opts.useTailApprox);
L{end+1} = '';
L{end+1} = '| Omnibus | Scope | real Q | p |';
L{end+1} = '|---|---|---|---|';
for o = stats.omnibus
    L{end+1} = sprintf('| %s | %s | %.1f | %.4g%s |', o.type, o.label, o.emp, o.p, local_tf(o.tailFitUsed)); %#ok<AGROW>
end
L{end+1} = '';
L{end+1} = '| Max statistic | Family | Direction | FWE threshold (R) | highest real R | Significant: single-step | step-down |';
L{end+1} = '|---|---|---|---|---|---|---|';
for e = stats.maxstat
    L{end+1} = sprintf('| %s | %s | %s | %.3f | %.3f | %d | %d |', e.label, e.family, e.tail, e.threshold, ...
        e.observedMax, e.NsigSingle, e.NsigStepdown); %#ok<AGROW>
end
if any(strcmp({stats.maxstat.family}, 'allanalyses')) && isfield(mapset.meta, 'settings')
    L{end+1} = '';
    L{end+1} = ['- ''allanalyses'' combines analyses that share their shuffles (subscores: same chains; ', ...
        'patients missing a score are skipped for that score, see shuffle jointness above).'];
end
ea_permtest_readme_append(mapset.readmeDir, 'Permutation statistics', L);


function s = local_tf(tf)
if tf
    s = ' (GPD)';
else
    s = '';
end


function opts = local_defaults(opts, defaults)
fn = fieldnames(defaults);
for i = 1:numel(fn)
    if ~isfield(opts, fn{i})
        opts.(fn{i}) = defaults.(fn{i});
    end
end


function s = local_esc(s)
% ea_error passes its message through sprintf -- keep paths/messages literal.
s = strrep(strrep(s, '\', '\\'), '%', '%%');
