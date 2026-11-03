function [permset, permsetfile] = ea_permtest_make_permset(groupdir, patientIDs, Nperm, rngseed)
% Create the cohort-level shuffle set ("permset") shared by all permutation
% tests in explorers/permutation_testing/, and save it to
% <groupdir>/permsets/.
%
% A permset does NOT store shuffled scores. For every permutation it stores
% a random ordering (rank) of ALL patients in the lead group.
% ea_permtest_derive_shuffles.m turns that into a shuffle of the scores
% among the patients actually selected for an analysis, by restricting the
% ordering to the selected patients. A random ordering of the whole cohort,
% restricted to any subset, is a uniformly random ordering of that subset,
% so one permset gives a valid permutation test for any patient selection,
% and identical selections always get identical shuffles (across
% modalities and explorer tools). Shuffling the whole cohort's SCORES and
% then subsetting would not be valid: selected patients would receive
% scores of patients that are not part of the analysis.
%
% groupdir   - lead group folder, i.e. fileparts(obj.leadgroup).
% patientIDs - Npatients x 1 cellstr of patient identifiers in lead group
%              (M.patient.list) order, as returned by the adapters.
% Nperm      - number of permutations. Default 1000.
% rngseed    - seed passed to rng(). Default 'default'.
%
% Extra orderings ("spares", 20% of Nperm, at least 100) are stored after
% the first Nperm. ea_permtest_derive_shuffles.m uses them, in order, to
% replace shuffles that turn out identical to an earlier one (or to the
% unshuffled data) once restricted to a small selection. Because columns
% are drawn sequentially from the same seed, the first Nperm orderings do
% not depend on the number of spares.
%
% Never overwrites an existing permset -- use ea_permtest_get_permset.m to
% load an existing one (or create it if missing).

if ~exist('Nperm', 'var') || isempty(Nperm)
    Nperm = 1000;
end
if ~exist('rngseed', 'var') || isempty(rngseed)
    rngseed = 'default';
end

if ~isscalar(Nperm) || Nperm < 1 || Nperm ~= round(Nperm)
    ea_error('Nperm must be a positive integer.');
end

patientIDs = cellstr(patientIDs(:));
Npatients = numel(patientIDs);
if Npatients < 2
    ea_error('A permset needs at least 2 patients.');
end
if numel(unique(patientIDs)) < Npatients
    ea_error('patientIDs must be unique -- cannot build a permset for a lead group with duplicate patient IDs.');
end

[fingerprint, filename] = ea_permtest_permset_fingerprint(patientIDs, Nperm, rngseed);
permsetdir = fullfile(groupdir, 'permsets');
permsetfile = fullfile(permsetdir, filename);
if exist(permsetfile, 'file')
    ea_error(['Permset already exists:\n', local_esc(permsetfile), '\n', ...
        'Use ea_permtest_get_permset to load it.']);
end

% Create the folder first (and fail with the actual reason, e.g. a lead
% group path that does not exist on this machine or a read-only drive).
if ~isfolder(groupdir)
    ea_error(['Lead group folder does not exist on this machine:\n', local_esc(groupdir), '\n', ...
        'Was the explorer loaded from another location, or is the drive not mounted?']);
end
if ~isfolder(permsetdir)
    [ok, msg] = mkdir(permsetdir);
    if ~ok
        ea_error(['Cannot create the permset folder:\n', local_esc(permsetdir), '\n', local_esc(msg)]);
    end
end

Nspare = max(100, ceil(0.2 * Nperm));
Ncols = Nperm + Nspare;

if Npatients <= intmax('uint16')
    rankClass = 'uint16';
else
    rankClass = 'uint32';
end

% Leave the session's global random stream as it was found.
prevState = rng;
rng(rngseed);
rngState = rng;
Ranks = zeros(Npatients, Ncols, rankClass);
for c = 1:Ncols
    Ranks(:, c) = randperm(Npatients)';
end
rng(prevState);

permset = struct;
permset.formatversion = 1;
permset.fingerprint = fingerprint;
permset.patientIDs = patientIDs;
permset.Npatients = Npatients;
permset.Nperm = Nperm;
permset.Nspare = Nspare;
permset.rngseed = rngseed;
permset.rngSeedUsed = rngState.Seed;
permset.rngType = rngState.Type;
permset.created = datestr(now, 'yyyy-mm-dd HH:MM:SS');
permset.matlabVersion = version;
permset.Ranks = Ranks; % Npatients x (Nperm+Nspare), column c = random rank of each patient

save(permsetfile, 'permset', '-v7');

fprintf('Created permset (%d patients, %d permutations + %d spares):\n%s\n', ...
    Npatients, Nperm, Nspare, permsetfile);


function s = local_esc(s)
% ea_error passes its message through sprintf -- keep paths/messages literal.
s = strrep(strrep(s, '\', '\\'), '%', '%%');
