function [Iperm, PermIdx, info] = ea_permtest_derive_shuffles(permset, I, sel, group)
% Derive the shuffles of the outcome(s) for one analysis from a
% cohort-level permset (see ea_permtest_make_permset.m).
%
% Each permset column is read as a chain of donations over the whole
% cohort: every patient receives the score of one donor (the patient
% ranked j-th receives the score of the j-th patient in lead group order).
% For an analysis, patients that do not take part -- not selected, or
% missing the score of that outcome -- are SKIPPED in the chain: whoever
% would have received from them receives from their donor instead (and so
% on, until a participating donor is reached). Skipping patients from a
% random chain leaves a random chain of the remaining patients, so every
% outcome gets a valid (uniformly random) shuffle of exactly its own
% patients, and participating patients keep their donor unless that donor
% was skipped.
%
% Example: chain A<-E, B<-B, C<-F, D<-A, E<-D, F<-C. Selection A,C,D,F
% (scores 10,20,30,40): A's donor E is skipped, so A receives from E's
% donor D -> A=30, C=40, D=10, F=20; B and E keep their own value.
%
% Several outcomes (e.g. subscores) are read from the SAME chain:
% - no missing scores: all outcomes get the identical shuffle, i.e. every
%   patient receives their whole row of scores from one donor;
% - a patient missing one outcome is skipped only for that outcome, so
%   that outcome keeps exactly its own set of patients (and E-fields) in
%   every permutation, and only the patient that would have received from
%   them gets a different donor for that outcome.
% Patients missing every outcome, and patients outside the selection,
% keep their own values. Shuffles are applied per patient, i.e. BEFORE any
% side expansion or mirroring, so mirrored copies stay tied.
%
% permset - struct from ea_permtest_get_permset / ea_permtest_make_permset.
% I       - Npatients x K outcome(s), one row per patient (lead group
%           order), one column per outcome.
% sel     - indices of selected patients (obj.patientselection). Order does
%           not matter. Default: all patients.
% group   - Npatients x 1 stratum label per patient (e.g.
%           obj.M.patient.group) to shuffle only within strata (each
%           stratum has its own, independent chain), or [] to shuffle
%           across all selected patients. Default [].
%
% Iperm   - Npatients x Nperm (K = 1) or Npatients x Nperm x K:
%           Iperm(:,p,k) == I(PermIdx(:,p,k), k).
% PermIdx - Npatients x Nperm (x K): the patient whose score patient i
%           receives in permutation p (for outcome k).
% info    - struct documenting how the shuffles were obtained:
%   .mode               'derived'   - random chains from the permset
%                       'exact'     - all possible arrangements enumerated
%                                     (selection too small for sampling)
%                       'enumerated-sampled' - Nperm distinct arrangements
%                                     sampled from the full enumeration
%   .NpermRequested     permset.Nperm
%   .Nperm              number of shuffles returned (< requested only for 'exact')
%   .Nreplaced          permset columns skipped because, for this selection,
%                       they duplicated an earlier shuffle or left every
%                       outcome unshuffled ('derived' only)
%   .columnsUsed        permset columns used, in order ('derived' only)
%   .totalArrangements  number of possible arrangements of the participating
%                       patients (all outcomes together)
%   .stratified, .strata, .selection, .Noutcomes
%   .poolSizes          Noutcomes x Nstrata, patients shuffled per outcome
%   .NexcludedNaN       1 x Noutcomes, selected patients missing that outcome
%   .NunshuffledPerOutcome / .NduplicatesPerOutcome
%                       1 x Noutcomes; nonzero only when outcomes have
%                       different missing patients and small pools
%   .jointness          mean fraction of patients with all outcomes whose
%                       scores all come from the same donor (1 = fully joint)
%   .permsetFingerprint
%
% The result depends only on (permset, NaN pattern of I, sel, group): the
% same selection always gets the same shuffles. In 'derived' mode the
% shuffles are distinct as a whole (all outcomes together) and never leave
% every outcome unshuffled. In 'exact' / 'enumerated-sampled' mode every
% chain over the participating patients is used at most once (the
% unshuffled one never), which keeps each outcome's shuffles exactly
% uniform even when outcomes have different missing patients.

Npatients = permset.Npatients;
if isvector(I)
    I = I(:);
end
if size(I, 1) ~= Npatients
    ea_error(sprintf('Outcome has %d rows but the permset covers %d patients.', size(I, 1), Npatients));
end
K = size(I, 2);

if ~exist('sel', 'var') || isempty(sel)
    sel = 1:Npatients;
end
if any(sel(:) < 1 | sel(:) > Npatients | sel(:) ~= round(sel(:)))
    ea_error('Patient selection contains invalid indices.');
end
sel = unique(sel(:))'; % sorted -- the order a selection is given in must not change the shuffles

stratified = exist('group', 'var') && ~isempty(group);
if stratified
    group = group(:);
    if numel(group) ~= Npatients
        ea_error(sprintf('group has %d entries but the permset covers %d patients.', numel(group), Npatients));
    end
    if any(isnan(group(sel)))
        ea_error('group labels of selected patients must not be NaN.');
    end
    strata = unique(group(sel))';
    universes = arrayfun(@(s) find(group == s)', strata, 'Uni', 0); % whole cohort's members of each stratum
else
    strata = 1;
    universes = {1:Npatients};
end
G = numel(strata);

% Participating patients per stratum and outcome
pools = cell(G, K);
for g = 1:G
    inStratum = ismember(sel, universes{g});
    for k = 1:K
        pools{g, k} = sel(inStratum & ~isnan(I(sel, k))');
    end
end
poolSizes = reshape(cellfun(@numel, pools), G, K)'; % K x G
unionPools = cell(1, G);
for g = 1:G
    unionPools{g} = unique([pools{g, :}]);
end
movable = cellfun(@numel, unionPools) > 1;
if ~any(movable)
    ea_error('Nothing to shuffle: no stratum contains at least 2 selected patients with a non-NaN score.');
end
shuffledOutcome = any(poolSizes > 1, 2)'; % outcomes with anything to shuffle
unionRows = [unionPools{movable}];

totalArrangements = prod(factorial(cellfun(@numel, unionPools(movable)))); % Inf if astronomically many
Nperm = permset.Nperm;

% Below this many possible arrangements, random sampling would mostly
% produce duplicates -- enumerate all arrangements instead.
maxEnumerate = 10 * Nperm;

Nreplaced = 0;
columnsUsed = [];
if totalArrangements <= maxEnumerate
    % All chains over the participating patients (the unshuffled one
    % removed). Every outcome's induced shuffles are then exactly uniform.
    sigma = local_enumerate(unionPools(movable), Npatients);
    sigma(:, all(sigma == (1:Npatients)', 1)) = [];
    PermAll = local_restrict_all(sigma, pools, K, Npatients);
    Ncand = size(sigma, 2);
    if Ncand <= Nperm
        PermIdx = PermAll;
        mode = 'exact';
        ea_cprintf('CmdWinWarnings', ['Only %d distinct shuffles exist for this selection; using all of them ', ...
            '(exact permutation test) instead of %d random permutations.\n'], Ncand, Nperm);
    else
        stream = RandStream('mt19937ar', 'Seed', permset.rngSeedUsed); % deterministic, leaves the global stream untouched
        pick = sort(randperm(stream, Ncand, Nperm));
        PermIdx = PermAll(:, pick, :);
        mode = 'enumerated-sampled';
        fprintf(['Only %d distinct shuffles exist for this selection; sampled %d distinct ones ', ...
            'from the full enumeration.\n'], Ncand, Nperm);
    end
else
    % Chains from the permset: within each stratum's universe, the patient
    % ranked j-th receives from the j-th member (lead group order).
    Ncols = size(permset.Ranks, 2);
    sigma = repmat((1:Npatients)', 1, Ncols);
    for g = 1:G
        U = universes{g}(:);
        [~, ord] = sort(permset.Ranks(U, :), 1);
        receivers = U(ord);
        sigma(receivers + (0:Ncols-1) * Npatients) = repmat(U, 1, Ncols);
    end
    PermAll = local_restrict_all(sigma, pools, K, Npatients);

    % Keep each distinct shuffle (all outcomes together) once, in permset
    % column order, dropping any that leave every outcome unshuffled
    % (the unshuffled arrangement is prepended as row 1).
    keys = reshape(permute(PermAll(unionRows, :, shuffledOutcome), [2 1 3]), Ncols, []);
    idKey = repmat(unionRows, 1, sum(shuffledOutcome));
    [~, first] = unique([idKey; keys], 'rows', 'stable');
    first = first(first > 1) - 1;
    if numel(first) < Nperm
        ea_error(sprintf(['Ran out of permset orderings: only %d distinct shuffles found for this ', ...
            'selection (%d needed). Use a smaller Nperm.'], numel(first), Nperm));
    end
    columnsUsed = first(1:Nperm)';
    PermIdx = PermAll(:, columnsUsed, :);
    Nreplaced = columnsUsed(end) - Nperm;
    mode = 'derived';
    if Nreplaced > 0
        fprintf(['%d permset ordering(s) gave a duplicate (or unshuffled) arrangement for this ', ...
            'selection and were replaced by spares.\n'], Nreplaced);
    end
end
NpermOut = size(PermIdx, 2);

Iperm = zeros(Npatients, NpermOut, K);
for k = 1:K
    Ik = I(:, k);
    Iperm(:, :, k) = Ik(PermIdx(:, :, k));
end

% Bookkeeping per outcome
NunshuffledPerOutcome = zeros(1, K);
NduplicatesPerOutcome = zeros(1, K);
for k = find(shuffledOutcome)
    rows = [pools{:, k}];
    Pk = PermIdx(rows, :, k)';
    NunshuffledPerOutcome(k) = sum(all(Pk == rows, 2));
    NduplicatesPerOutcome(k) = NpermOut - size(unique(Pk, 'rows'), 1);
end
complete = sel(all(~isnan(I(sel, :)), 2)');
if K > 1 && ~isempty(complete)
    sameDonor = all(PermIdx(complete, :, :) == PermIdx(complete, :, 1), 3);
    jointness = mean(sameDonor(:));
else
    jointness = 1;
end

if K == 1
    PermIdx = PermIdx(:, :, 1);
    Iperm = Iperm(:, :, 1);
end

info = struct;
info.mode = mode;
info.NpermRequested = Nperm;
info.Nperm = NpermOut;
info.Nreplaced = Nreplaced;
info.columnsUsed = columnsUsed;
info.totalArrangements = totalArrangements;
info.stratified = stratified;
info.strata = strata;
info.poolSizes = poolSizes;
info.selection = sel;
info.Noutcomes = K;
info.NexcludedNaN = numel(sel) - sum(poolSizes, 2)';
info.NunshuffledPerOutcome = NunshuffledPerOutcome;
info.NduplicatesPerOutcome = NduplicatesPerOutcome;
info.jointness = jointness;
info.permsetFingerprint = permset.fingerprint;


function PermAll = local_restrict_all(sigma, pools, K, Npatients)
% Restrict the chains (receiver -> donor, Npatients x Ncols) to each
% outcome's participating patients by skipping non-participants.
Ncols = size(sigma, 2);
PermAll = repmat((1:Npatients)', 1, Ncols, K);
for k = 1:K
    for g = 1:size(pools, 1)
        P = pools{g, k};
        if numel(P) > 1
            PermAll(P, :, k) = local_first_return(sigma, P(:), Npatients);
        end
    end
end


function D = local_first_return(sigma, P, Npatients)
% For each participant i in P (and chain column c): follow the chain
% i <- sigma(i) <- sigma(sigma(i)) ... until the first participant donor.
% Terminates because every chain is a cycle through i itself.
inP = false(Npatients, 1);
inP(P) = true;
Ncols = size(sigma, 2);
D = sigma(P, :);
colIdx = repmat(1:Ncols, numel(P), 1);
skip = ~inP(D);
while any(skip(:))
    D(skip) = sigma(D(skip) + (colIdx(skip) - 1) * Npatients);
    skip = ~inP(D);
end


function sigma = local_enumerate(pools, Npatients)
% All possible chains over the given pools (cartesian product across
% strata), including the unshuffled one. Column m is receiver -> donor.
G = numel(pools);
P = cell(1, G);
nArr = zeros(1, G);
for g = 1:G
    P{g} = perms(pools{g}); % rows = all orderings of the pool's patients
    nArr(g) = size(P{g}, 1);
end

if G == 1
    grids = {(1:nArr(1))'};
else
    ranges = arrayfun(@(n) 1:n, nArr, 'Uni', 0);
    grids = cell(1, G);
    [grids{:}] = ndgrid(ranges{:});
end

sigma = repmat((1:Npatients)', 1, prod(nArr));
for g = 1:G
    sigma(pools{g}, :) = P{g}(grids{g}(:), :)';
end
