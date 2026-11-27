function [permset, permsetfile] = ea_permtest_get_permset(groupdir, patientIDs, Nperm, rngseed, createIfMissing)
% Load the cohort-level shuffle set matching (patientIDs, Nperm, rngseed)
% from <groupdir>/permsets/, creating it if it does not exist yet. See
% ea_permtest_make_permset.m for what a permset contains.
%
% The loaded file is verified against the request (fingerprint, patient
% list and order, Nperm, size of the ordering matrix), so a permset is
% never silently used for a different cohort.
%
% groupdir        - lead group folder, i.e. fileparts(obj.leadgroup).
% patientIDs      - Npatients x 1 cellstr in lead group order.
% Nperm           - number of permutations. Default 1000.
% rngseed         - seed. Default 'default'.
% createIfMissing - create the permset if no matching file exists.
%                   Default true.

if ~exist('Nperm', 'var') || isempty(Nperm)
    Nperm = 1000;
end
if ~exist('rngseed', 'var') || isempty(rngseed)
    rngseed = 'default';
end
if ~exist('createIfMissing', 'var') || isempty(createIfMissing)
    createIfMissing = true;
end

patientIDs = cellstr(patientIDs(:));
[fingerprint, filename] = ea_permtest_permset_fingerprint(patientIDs, Nperm, rngseed);
permsetfile = fullfile(groupdir, 'permsets', filename);

if ~exist(permsetfile, 'file')
    if createIfMissing
        [permset, permsetfile] = ea_permtest_make_permset(groupdir, patientIDs, Nperm, rngseed);
        return
    else
        ea_error(sprintf('No permset found for this cohort/Nperm/seed:\n%s', permsetfile));
    end
end

S = load(permsetfile, 'permset');
if ~isfield(S, 'permset')
    ea_error(sprintf('File does not contain a permset:\n%s', permsetfile));
end
permset = S.permset;

if ~strcmp(permset.fingerprint, fingerprint) ...
        || ~isequal(permset.patientIDs, patientIDs) ...
        || permset.Nperm ~= Nperm ...
        || ~isequal(size(permset.Ranks), [numel(patientIDs), permset.Nperm + permset.Nspare])
    ea_error(sprintf(['Permset file does not match the requested cohort/Nperm/seed ', ...
        '(file changed or corrupted?):\n%s'], permsetfile));
end

fprintf('Using permset (%d patients, %d permutations):\n%s\n', ...
    permset.Npatients, permset.Nperm, permsetfile);
