function [pval, tailFitUsed, apar, kpar, upar] = ea_unified_perm_palm_gpd_pval(nullvals, empval, tail, Pthr)
% Permutation p-value with an optional Generalized Pareto Distribution
% (GPD) tail refinement, ported from FSL PALM's palm_pareto.m (Winkler et
% al. 2016, NeuroImage 141:502-516; method originates from Knijnenburg et
% al. 2009, Bioinformatics 25(12):i161-8; GPD fit from Hosking & Wallis
% 1987, Technometrics 29:339-349; goodness-of-fit test from Choulakian &
% Stephens 2001, Technometrics 43(4):478-484). Self-contained -- does not
% require the separate PALM toolbox on the MATLAB path.
%
% If you use this refinement (useTailApprox=true and the fit is actually
% applied) in a publication, cite Winkler et al. 2016 (PALM) and
% Knijnenburg et al. 2009 (the method's origin). A popup reminder is shown
% the first time a fit is accepted in a session (see below), falling back
% to a console message in headless/batch environments.
%
% Why: the plain rank-based permutation p-value (see
% ea_unified_perm_nulldist_plot) has a resolution floor of 1/(Nperm+1) -- it
% cannot express anything smaller no matter how extreme the empirical
% value truly is. When the rank-based p is already small (below Pthr),
% this fits a smooth GPD curve to the shape of the null distribution's
% extreme tail and reads the p-value off that curve instead, which can
% resolve p-values below the permutation-count floor. The fit is only used
% if it actually looks GPD-shaped (checked via an Anderson-Darling
% goodness-of-fit test at each candidate tail cutoff); otherwise this
% silently falls back to the plain rank-based p, so a bad fit never
% fabricates a number.
%
% Scope: this operates on ONE shared null distribution and a single
% empirical statistic (e.g. the max-stat or Eisenstein-Q omnibus null),
% matching how ea_unified_perm_nulldist_plot is used. It is deliberately not
% applied to the per-voxel uncorrected test, where every voxel has its own
% separate null and a GPD fit would have to be repeated per voxel --
% multiplying cost back up for a test where this level of precision isn't
% the point.
%
% nullvals - Nperm x 1 permutation null distribution.
% empval   - the empirical (observed, unpermuted) statistic.
% tail     - 'right' (large nullvals are extreme, e.g. positive max-stat/Eisenstein Q)
%            'left'  (small nullvals are extreme, e.g. negative max-stat)
%            'both'  (large |nullvals| are extreme, e.g. out-of-sample prediction R)
% Pthr     - only p-values below this are candidates for tail refinement
%            (matches PALM's -accel tail default of 0.10). Default 0.10.
%
% Returns:
% pval        - the (possibly tail-refined) p-value. Uses the same
%               (count+1)/(Nperm+1) convention as ea_unified_perm_nulldist_plot
%               when no tail fit is used or applicable -- never exactly 0.
% tailFitUsed - true if a GPD tail fit passed its goodness-of-fit check AND
%               the empirical value actually fell inside that fitted tail.
% apar, kpar, upar - GPD scale, shape, and threshold (location) parameters
%               of the accepted fit, or NaN if no fit was used.

if ~exist('Pthr', 'var') || isempty(Pthr)
    Pthr = 0.10;
end

nullvals = nullvals(:);
nP = numel(nullvals); % NaN permutations stay IN this count (see substitution below)

switch tail
    case 'right'
        rev = false;
        g = nullvals;
        e = empval;
    case 'left'
        rev = true;
        g = nullvals;
        e = empval;
    case 'both'
        rev = false;
        g = abs(nullvals);
        e = abs(empval);
    otherwise
        error('ea_unified_perm_palm_gpd_pval:tail', 'tail must be ''right'', ''left'', or ''both''.');
end

% A permutation that produced no defined statistic (NaN) is evidence
% against the null, not missing data -- same convention as
% ea_unified_perm_nulldist_plot. It must never be counted as "exceeding" and
% must never enter the fitted tail, but it does stay in the denominator
% (nP above). Substituting the least-extreme-possible value for the
% current tail direction achieves both at once: it sorts to the bulk end,
% so it can never be selected into Gtail below, while still occupying one
% slot of nP.
if rev
    g(isnan(g)) = +Inf; % 'left': small values are extreme, so NaN -> least extreme = +Inf
else
    g(isnan(g)) = -Inf; % 'right'/'both': large values are extreme, so NaN -> least extreme = -Inf
end

% Plain rank-based p-value first (same convention as ea_unified_perm_nulldist_plot).
if isnan(e) || nP == 0
    pval = NaN;
    tailFitUsed = false;
    apar = NaN; kpar = NaN; upar = NaN;
    return
end
if rev
    exceedCount = sum(g <= e);
else
    exceedCount = sum(g >= e);
end
pval = (exceedCount + 1) / (nP + 1);

tailFitUsed = false;
apar = NaN; kpar = NaN; upar = NaN;

if pval >= Pthr
    return
end

% ---- GPD tail refinement (ported from PALM's palm_pareto.m) ----
if rev
    [Gsorted, tailFrac] = ea_unified_perm_ecdf(g, 'descend');
else
    [Gsorted, tailFrac] = ea_unified_perm_ecdf(g, 'ascend');
end

Q = (751:10:999)/1000; % candidate tail cutoffs (top 24.9% down to top 0.1%), most inclusive first
nQ = numel(Q);
q = 1;
fitAccepted = false;
while ~fitAccepted && q < nQ - 1
    qidx = tailFrac >= Q(q);
    Gtail = Gsorted(qidx);
    qi = find(qidx, 1);
    if qi == 1
        thisUpar = Gsorted(qi) - mean(Gsorted(qi:min(qi+1,nP)));
    else
        thisUpar = mean(Gsorted(qi-1:qi));
    end
    if rev
        ytail = thisUpar - Gtail;
    else
        ytail = Gtail - thisUpar;
    end

    % Method-of-moments GPD fit, Hosking & Wallis (1987) sect. 3.2.
    x = mean(ytail);
    s2 = var(ytail);
    if s2 == 0 || isnan(s2)
        q = q + 1;
        continue
    end
    thisApar = x*(x^2/s2 + 1)/2;
    thisKpar = (x^2/s2 - 1)/2;

    A2pval = ea_unified_perm_andersondarling(ea_unified_perm_gpdpval(ytail, thisApar, thisKpar), thisKpar);

    if A2pval > 0.05
        fitAccepted = true;
        apar = thisApar; kpar = thisKpar; upar = thisUpar;
        if rev
            qualifies = e < upar;
            y = upar - e;
        else
            qualifies = e > upar;
            y = e - upar;
        end
        if qualifies
            cte = numel(Gtail) / nP;
            pval = cte * ea_unified_perm_gpdpval(y, apar, kpar);
            tailFitUsed = true;
            ea_unified_perm_gpd_cite();
        end
        % if ~qualifies: a good GPD fit was found, but this particular
        % empirical value doesn't fall inside the fitted tail -- keep the
        % plain rank-based p-value computed above.
    else
        q = q + 1;
    end
end
end

% ==============================================================
function ea_unified_perm_gpd_cite()
% Shows a citation reminder popup once per MATLAB session, only when the
% GPD tail fit is actually applied to a reported p-value (not merely
% requested-but-inapplicable). Falls back to a console message in
% headless/batch environments (e.g. running on a cluster with no figure
% window support), where a popup would either error or block silently.
persistent shown
if ~isempty(shown)
    return
end
shown = true;

msg = sprintf(['A GPD tail-fit p-value was used in this analysis.\n\n', ...
    'If reporting this, please cite:\n\n', ...
    '1. Winkler AM, Ridgway GR, Douaud G, Nichols TE, Smith SM (2016).\n', ...
    '   Faster permutation inference in brain imaging.\n', ...
    '   NeuroImage 141:502-516. (FSL PALM)\n\n', ...
    '2. Knijnenburg TA, Wessels LFA, Reinders MJT, Shmulevich I (2009).\n', ...
    '   Fewer permutations, more accurate P-values.\n', ...
    '   Bioinformatics 25(12):i161-8.']);

if usejava('jvm') && feature('ShowFigureWindows')
    msgbox(msg, 'Citation reminder: GPD tail-fit p-value', 'help');
else
    fprintf('\n[ea_unified_perm_palm_gpd_pval] %s\n\n', strrep(msg, newline, sprintf('\n  ')));
end
end

% ==============================================================
function [Gsorted, tailFrac] = ea_unified_perm_ecdf(g, ord)
% Empirical CDF with modified-competition-ranking tie handling (tied
% values all share the most-inclusive count), matching PALM's
% palm_competitive(...,true) as used inside palm_pareto.m.
nP = numel(g);
Gsorted = sort(g, ord);
tailFrac = zeros(nP, 1);
i = 1;
while i <= nP
    j = i;
    while j < nP && Gsorted(j+1) == Gsorted(i)
        j = j + 1;
    end
    tailFrac(i:j) = j / nP;
    i = j + 1;
end
end

% ==============================================================
function p = ea_unified_perm_gpdpval(x, a, k)
% GPD survival function (Hosking & Wallis 1987 parameterization).
if abs(k) < eps
    p = exp(-x/a);
else
    p = (1 - k*x/a).^(1/k);
end
if k > 0
    p(x > a/k) = 0;
end
end

% ==============================================================
function A2pval = ea_unified_perm_andersondarling(z, k)
% Anderson-Darling goodness-of-fit p-value for a GPD fit, using the
% interpolation table from Choulakian & Stephens (2001), Table 2 (Case 3:
% both scale and shape estimated). Ported directly from PALM's palm_pareto.m.
ktable = [0.9 0.5 0.2 0.1 0 -0.1 -0.2 -0.3 -0.4 -0.5]';
ptable = [0.5 0.25 0.1 0.05 0.025 0.01 0.005 0.001];
A2table = [ ...
    0.3390 0.4710 0.6410 0.7710 0.9050 1.0860 1.2260 1.5590
    0.3560 0.4990 0.6850 0.8300 0.9780 1.1800 1.3360 1.7070
    0.3760 0.5340 0.7410 0.9030 1.0690 1.2960 1.4710 1.8930
    0.3860 0.5500 0.7660 0.9350 1.1100 1.3480 1.5320 1.9660
    0.3970 0.5690 0.7960 0.9740 1.1580 1.4090 1.6030 2.0640
    0.4100 0.5910 0.8310 1.0200 1.2150 1.4810 1.6870 2.1760
    0.4260 0.6170 0.8730 1.0740 1.2830 1.5670 1.7880 2.3140
    0.4450 0.6490 0.9240 1.1400 1.3650 1.6720 1.9090 2.4750
    0.4680 0.6880 0.9850 1.2210 1.4650 1.7990 2.0580 2.6740
    0.4960 0.7350 1.0610 1.3210 1.5900 1.9580 2.2430 2.9220];

k = max(0.5, k);
z = flipud(z(:))';
n = numel(z);
j = 1:n;

A2 = -n - (1/n)*((2*j-1)*(log(z) + log(1 - z(n+1-j)))');
i1 = interp1(ktable, A2table, k, 'linear', 'extrap');
i2 = interp1(i1, ptable, A2, 'linear', 'extrap');
A2pval = max(min(i2, 1), 0);
end
