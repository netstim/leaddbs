function ea_unified_perm_sweetspot_consistency_check(obj, patsel, Remp)
% Cross-checks ea_unified_perm_sweetspot_calcstats.m's empirical (real,
% unpermuted) map against the live, shared ea_unifiedmapping_calcstats.m
% computing the same thing -- both implement the same sweetspotmapping +
% Correlations coverage-masking logic, deliberately kept separate (see
% ea_unified_perm_sweetspot_calcstats.m's header for why), which means they
% can drift apart if the shared function is ever changed without updating
% this feature's copy. This check catches that drift at the start of every
% permtest run, before the expensive Nperm permutation loop, rather than
% silently producing a subtly-wrong empirical/null map.
%
% obj/patsel - same object/patient selection permtest was called with.
% Remp       - the {side} cell of empirical R-values already computed by
%              ea_unified_perm_sweetspot_calcstats (skipsigthresh=true,
%              i.e. unmasked).
%
% Does not error or block execution -- mismatches are surfaced as a loud
% (popup, with console fallback) warning, since permtest's own Remp is
% still usable even if this check fails; this is a correctness alarm for
% you to act on, not a hard gate.

% Live calcstats applies significance masking when obj.showsignificantonly
% is on; Remp is deliberately unmasked, so mask must be off for a fair
% comparison here. onCleanup guarantees restoration even if the call below
% errors.
origSig = obj.showsignificantonly;
restoreSig = onCleanup(@() setshowsig(obj, origSig)); %#ok<NASGU>
obj.showsignificantonly = false;

liveVals = ea_unifiedmapping_calcstats(obj, patsel);

tol = 1e-6;
mismatchReport = {};
for side = 1:numel(Remp)
    a = Remp{side}(:);
    b = liveVals{side}(:);
    if numel(a) ~= numel(b)
        mismatchReport{end+1} = sprintf('side %d: voxel count differs (%d vs %d)', side, numel(a), numel(b)); %#ok<AGROW>
        continue
    end
    mismatch = (abs(a-b) > tol) | (isnan(a) ~= isnan(b));
    nMismatch = sum(mismatch);
    if nMismatch > 0
        mismatchReport{end+1} = sprintf('side %d: %d/%d voxels differ (max |diff| = %.4g)', ...
            side, nMismatch, numel(a), max(abs(a(~isnan(a)&~isnan(b)) - b(~isnan(a)&~isnan(b))), [], 'omitnan')); %#ok<AGROW>
    end
end

if isempty(mismatchReport)
    return
end

msg = sprintf(['ea_unified_perm_sweetspot_calcstats.m''s empirical map does NOT match\n', ...
    'ea_unifiedmapping_calcstats.m computing the same (real, unpermuted) data.\n\n', ...
    'This means the two implementations have drifted apart -- most likely\n', ...
    'ea_unifiedmapping_calcstats.m was changed without mirroring the change into\n', ...
    'ea_unified_perm_sweetspot_calcstats.m. permtest results below may be wrong.\n\n', ...
    '%s'], strjoin(mismatchReport, sprintf('\n')));

warning('ea_unified_perm_sweetspot:calcstats_drift', '%s', msg);
if usejava('jvm') && feature('ShowFigureWindows')
    msgbox(msg, 'permtest consistency check failed', 'warn');
end
end

function setshowsig(obj, val)
obj.showsignificantonly = val;
end
