function [run, bundle, report] = ea_permtest_fibers_disctract(obj, opts)
% One-call fiber permutation test for the fiber filtering explorer
% (ea_disctract): adapter -> consistency check (stops if it fails) ->
% cohort shuffle set -> shuffles -> permuted fiber maps saved to a new run
% folder (see ea_permtest_fibers_run.m for the folder contents).
%
% obj  - ea_disctract object (e.g. loaded from the .fibfilt; obj is left
%        as it was found).
% opts - struct, all optional:
%   .fibfiltPath      the .fibfilt obj was loaded from; the lead group
%                     folder is taken from its location (recommended when
%                     working from the command line, as paths stored in
%                     the explorer can be outdated)
%   .groupdir         lead group folder directly (overrides fibfiltPath)
%   .name             name of the analysis: results go to
%                     <lead group>/fiberfiltering/permutation/<name>/fibers/
%                     (see ea_permtest_setup.m). Default: analysis_<time>.
%   .Nperm            number of permutations. Default 1000.
%   .stratifyByGroup  shuffle only within cohort groups (M.patient.group).
%                     Default false.
%   .fastPath         fast computation for Spearman / Pearson (checked
%                     against the test function). Default true.
%   .maxWorkers       'auto' (default: cores - 1, at most 4) or a number;
%                     0/1 runs serially. Only used without the fast path.
%                     Results do not depend on it.
%   .batchSize        permutations per write to disk. Default 100.
%
% run    - see ea_permtest_fibers_run.m (run.rundir is the run folder).
% bundle - the adapter bundle used.
% report - the consistency-check report.
%
% Example:
%   S = load(fibfiltPath, '-mat'); obj = S.tractset;
%   run = ea_permtest_fibers_disctract(obj, struct('fibfiltPath', fibfiltPath, 'name', 'OCD_N32_pooled'));

if ~exist('opts', 'var') || isempty(opts)
    opts = struct;
end
pathopts = struct;
for f = {'fibfiltPath', 'groupdir'}
    if isfield(opts, f{1})
        pathopts.(f{1}) = opts.(f{1});
    end
end
Nperm = 1000;
if isfield(opts, 'Nperm') && ~isempty(opts.Nperm)
    Nperm = opts.Nperm;
end

fprintf('Collecting fiber data from the explorer...\n');
bundle = ea_permtest_adapter_disctract(obj, pathopts);

fprintf('Checking that the permutation engine reproduces the explorer...\n');
report = ea_permtest_check_disctract(obj, bundle);
if ~report.ok
    ea_error(['The permutation engine does not reproduce the explorer''s fiber values for this analysis ', ...
        '(see the check output above). No permutations were calculated.'], 'title', 'Permutation testing');
end

[permset, permsetFile] = ea_permtest_get_permset(bundle.groupdir, bundle.patientIDs, Nperm, bundle.rngseed);

runopts = struct('checkReport', report, 'permsetFile', permsetFile);
for f = {'name', 'stratifyByGroup', 'fastPath', 'maxWorkers', 'batchSize'}
    if isfield(opts, f{1})
        runopts.(f{1}) = opts.(f{1});
    end
end
run = ea_permtest_fibers_run(bundle, permset, runopts);
