function [fingerprint, filename, seedstr] = ea_permtest_permset_fingerprint(patientIDs, Nperm, rngseed)
% Identity of a cohort-level shuffle set (permset): which patients it
% covers (and in which order), how many permutations, and which seed.
% Two calls with the same inputs always give the same fingerprint and file
% name, so any tool/modality asking for the same permset finds the same
% file in <leadgroup dir>/permsets/.
%
% patientIDs - Npatients x 1 cellstr, in lead group (M.patient.list) order.
% Nperm      - number of permutations.
% rngseed    - seed as passed to rng() (numeric or char, e.g. 'default').
%
% fingerprint - MD5 hex string.
% filename    - permset_N<Nperm>_seed<seed>_<first 8 chars of fingerprint>.mat
% seedstr     - rngseed as char.

if isnumeric(rngseed)
    seedstr = sprintf('%d', rngseed);
elseif ischar(rngseed) || isstring(rngseed)
    seedstr = char(rngseed);
else
    ea_error('rngseed must be numeric or char.');
end

key = struct('formatversion', 1, ...
    'patientIDs', {cellstr(patientIDs(:))}, ...
    'Nperm', double(Nperm), ...
    'rngseed', seedstr);

fingerprint = ea_datahash(key);
filename = sprintf('permset_N%d_seed%s_%s.mat', Nperm, regexprep(seedstr, '[^\w-]', ''), fingerprint(1:8));
