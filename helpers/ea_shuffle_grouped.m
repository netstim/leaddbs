function [Iperm, PermIdx] = ea_shuffle_grouped(I, Nperm, Sel, group, rngseed)
% Generate N random permutations of column vector I (one entry per patient),
% restricted to shuffle only within group(Sel) buckets (group-restricted
% permutation, e.g. do not mix scores between two cohorts/groups).
%
% Patients not in Sel are left untouched. Patients in Sel with a NaN score
% are excluded from the shuffle pool (kept in place); they are never used
% as a donor for another patient's slot either.
%
% I        - Npatients x 1 response variable (one score per patient, not yet
%            expanded to sides/mirrors). Hemiscore (2-column) I is not supported.
% Nperm    - number of permutations to generate
% Sel      - indices of patients eligible to participate in the shuffle
%            (typically obj.patientselection). Defaults to all patients.
% group    - Npatients x 1 group label per patient (e.g. obj.M.patient.group).
%            Pass [] to disable group restriction (all of Sel is one group).
% rngseed  - rng seed, defaults to 'default'.
%
% Returns:
% Iperm    - Npatients x Nperm, column p is I permuted for permutation p.
% PermIdx  - Npatients x Nperm, the patient-index mapping used for each
%            permutation (Iperm(:,p) == I(PermIdx(:,p))), saved for
%            reproducibility/auditing.
%
% Each new shuffle is checked against every previously accepted one and
% redrawn on a match, so all Nperm columns of PermIdx are guaranteed
% distinct (matches PALM's default Monte Carlo permutation behavior,
% which does the same unless run with -cmcp). In practice this rarely
% changes anything: with realistic patient counts, the total number of
% possible shuffles vastly exceeds any typical Nperm, so collisions are
% already vanishingly unlikely -- this just makes that guarantee exact
% rather than probabilistic, and only costs a redraw in the rare case it's
% needed.

if size(I,2) > 1
    error('ea_shuffle_grouped:hemiscore', 'Hemiscore responsevar (2 columns) is not yet supported for permutation testing.');
end

Npatients = length(I);

if ~exist('Sel','var') || isempty(Sel)
    Sel = 1:Npatients;
end
Sel = Sel(:)';

if ~exist('group','var') || isempty(group)
    group = ones(Npatients,1);
end
group = group(:);

if ~exist('rngseed','var') || isempty(rngseed)
    rngseed = 'default';
end

rng(rngseed);

PermIdx = repmat((1:Npatients)', 1, Nperm);
Iperm = repmat(I, 1, Nperm);

groups = unique(group(Sel));

maxAttempts = 1000; % redraw-on-duplicate cap; only matters if the shuffle pool is small relative to Nperm
warnedDuplicateLimit = false;

for p = 1:Nperm
    attempt = 0;
    while true
        idx = (1:Npatients)';
        for g = groups'
            pool = Sel(group(Sel)==g & ~isnan(I(Sel)));
            if numel(pool) > 1
                idx(pool) = pool(randperm(numel(pool)));
            end
        end
        attempt = attempt + 1;
        if p == 1 || ~any(all(idx == PermIdx(:,1:p-1), 1))
            break % unique among all previously accepted shuffles (or nothing to compare against yet)
        end
        if attempt >= maxAttempts
            if ~warnedDuplicateLimit
                warning('ea_shuffle_grouped:duplicateLimit', ...
                    ['Could not find a new, not-yet-used shuffle within %d attempts ', ...
                    '(permutation %d/%d) -- the shuffle pool may be small relative to Nperm. ', ...
                    'Accepting a duplicate here and suppressing this warning for the rest of this call.'], ...
                    maxAttempts, p, Nperm);
                warnedDuplicateLimit = true;
            end
            break
        end
    end
    PermIdx(:,p) = idx;
    Iperm(:,p) = I(idx);
end
