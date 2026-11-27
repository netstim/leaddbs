function shuffles = ea_permtest_bundle_shuffles(permset, bundle, stratifyByGroup)
% Shuffles for every analysis of an adapter bundle, derived from the
% cohort-level permset (see ea_permtest_derive_shuffles.m).
%
% Analyses that share the same patient selection (e.g. the subscores of
% 'Split & Color By Subscore') get ONE joint derivation over all their
% outcome columns, read from the same chains: without missing scores a
% patient's whole row of subscores moves together in every permutation;
% a patient missing a subscore is skipped for that subscore only (see
% ea_permtest_derive_shuffles.m). Analyses on different patients (e.g.
% 'Split & Color By Group') each get their own derivation, so a group's
% shuffles are identical to running that group's patients as a single
% analysis.
%
% permset         - from ea_permtest_get_permset.
% bundle          - from an adapter (ea_permtest_adapter_*.m).
% stratifyByGroup - shuffle only within bundle.group strata. Default false.
%
% shuffles(a), one per bundle.analyses(a):
%   .PermIdx    Npatients x Nperm, patient whose scores patient i receives
%   .Iperm      Npatients x Nperm, the analysis' outcome column shuffled
%   .derivation index of the joint derivation this analysis belongs to
%   .outcomePos position of this analysis' outcome within that derivation
%               (index into per-outcome info fields, e.g. info.NexcludedNaN)
%   .info       info struct of that derivation (see derive_shuffles)

if ~exist('stratifyByGroup', 'var') || isempty(stratifyByGroup)
    stratifyByGroup = false;
end
if stratifyByGroup
    strata = bundle.group;
else
    strata = [];
end

Nanalyses = numel(bundle.analyses);
selKeys = arrayfun(@(an) mat2str(sort(an.patientselection)), bundle.analyses, 'Uni', 0);
[~, ~, derivationOf] = unique(selKeys, 'stable');

shuffles = struct('PermIdx', cell(1, Nanalyses), 'Iperm', [], 'derivation', [], 'outcomePos', [], 'info', []);
for d = 1:max(derivationOf)
    members = find(derivationOf == d);
    cols = unique([bundle.analyses(members).outcomeCol], 'stable');
    sel = bundle.analyses(members(1)).patientselection;
    [~, PermIdx, info] = ea_permtest_derive_shuffles(permset, bundle.outcomes(:, cols), sel, strata);
    for a = members(:)'
        I = bundle.outcomes(:, bundle.analyses(a).outcomeCol);
        pos = find(cols == bundle.analyses(a).outcomeCol, 1);
        PermIdxA = PermIdx(:, :, pos);
        shuffles(a).PermIdx = PermIdxA;
        shuffles(a).Iperm = I(PermIdxA);
        shuffles(a).derivation = d;
        shuffles(a).outcomePos = pos;
        shuffles(a).info = info;
    end
end
