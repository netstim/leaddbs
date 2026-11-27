function outcomein = ea_permtest_fibers_outcome(I, bundle, a)
% Turn a per-patient outcome vector (original or one shuffled column) into
% the per-side outcome vectors that are fed into the fiber statistic of
% analysis a, exactly as the explorers' calcstats do it: duplicate to both
% sides, stack mirrored copies below the originals, then pick the
% analysis' selected (and mirrored) patients.
%
% Covariates are not supported yet: the adapters refuse to build a bundle
% when covariates are selected, and this function refuses a bundle that
% carries any, so covariate handling can never be applied half-way.
%
% I      - Npatients x 1 outcome (one column of bundle.outcomes, possibly
%          shuffled).
% bundle - input bundle from an adapter (ea_permtest_adapter_*.m).
% a      - index into bundle.analyses (subscore / group / single analysis).
%
% outcomein - 1 x Nsides cell, outcomein{side} is
%             numel(bundle.analyses(a).gpatsel) x 1.

if isfield(bundle, 'covars') && ~isempty(bundle.covars)
    ea_error(['Permutation testing is not designed for covariates yet. ', ...
        'No permutations were calculated.'], 'title', 'Permutation testing');
end
if ~isvector(I)
    ea_error('ea_permtest_fibers_outcome expects one outcome column at a time.');
end
I = I(:);
I = [I, I]; % both sides the same

if bundle.mirrorsides
    I = [I; I];
end

gpatsel = bundle.analyses(a).gpatsel;
outcomein = cell(1, bundle.Nsides);
for side = 1:bundle.Nsides
    outcomein{side} = I(gpatsel, side);
end
