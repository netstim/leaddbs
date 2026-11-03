function [run, cvb, report] = ea_permtest_cv_disctract(obj, opts)
% One-call prediction permutation test for the fiber filtering explorer
% (ea_disctract): adapter -> check against the explorer's own k-fold CV
% (stops if it fails) -> cohort shuffle set -> k-fold CV for every shuffle
% (see ea_permtest_cv_run.m for what is tested and saved).
%
% obj  - ea_disctract object (left as it was found).
% opts - struct, all optional:
%   .fibfiltPath / .groupdir  lead group folder, see ea_permtest_adapter_disctract.m
%   .name             analysis name: results go to
%                     <lead group>/fiberfiltering/permutation/<name>/prediction_k<k>/
%                     (an existing analysis with the same settings, e.g.
%                     with a fiber run, gets the prediction run added)
%   .kIter            k-fold iterations. Default obj.kIter.
%   .posthoccorrect   'Post-hoc correct for group' of the explorer. Default false.
%   .targetLabel      clinical variable for the combined prediction of all
%                     subscore models (optional, subscore mode only).
%   .Nperm            number of permutations. Default 1000.
%   .stratifyByGroup  shuffle only within cohort groups. Default false.
%   .confirmFcn       @(message) -> true/false, asked before long runs
%                     (k-fold iterations > 1).
%   .batchSize        shuffles computed together. Default 50.
%
% Example:
%   S = load(fibfiltPath, '-mat'); obj = S.tractset;
%   run = ea_permtest_cv_disctract(obj, struct('fibfiltPath', fibfiltPath, 'name', 'OCD_N32_pooled'));

if ~exist('opts', 'var') || isempty(opts)
    opts = struct;
end
pathopts = local_pick(opts, {'fibfiltPath', 'groupdir'});
cvopts = local_pick(opts, {'kIter', 'posthoccorrect', 'targetLabel'});
Nperm = 1000;
if isfield(opts, 'Nperm') && ~isempty(opts.Nperm)
    Nperm = opts.Nperm;
end

fprintf('Collecting the cross-validation setup from the explorer...\n');
cvb = ea_permtest_cv_adapter_disctract(obj, pathopts, cvopts);
report = ea_permtest_cv_check_disctract(obj, cvb);
if ~report.ok
    ea_error(['The permutation engine does not reproduce the explorer''s k-fold cross-validation ', ...
        'for these settings (see the command window). No permutations were calculated.'], 'title', 'Permutation testing');
end

[permset, permsetFile] = ea_permtest_get_permset(cvb.groupdir, cvb.patientIDs, Nperm, cvb.rngseed);
runopts = local_pick(opts, {'name', 'stratifyByGroup', 'confirmFcn', 'batchSize'});
runopts.checkReport = report;
runopts.permsetFile = permsetFile;
run = ea_permtest_cv_run(cvb, permset, runopts);


function s = local_pick(opts, fields)
s = struct;
for f = fields
    if isfield(opts, f{1})
        s.(f{1}) = opts.(f{1});
    end
end
