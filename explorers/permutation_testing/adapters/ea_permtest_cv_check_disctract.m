function report = ea_permtest_cv_check_disctract(obj, cvb, tol)
% Verify that ea_permtest_cv_predict.m reproduces the fiber filtering
% explorer's own k-fold cross-validation on the real scores: runs
% obj.kfoldcv (the code behind the explorer's 'k-fold (randomized)' CV)
% and compares its predictions patient by patient with the engine's.
% Run before any prediction permutation test.
%
% The explorer's CV is run as the explorer runs it ('Run CV' with the
% k-fold strategy: no custom selection, no imported model, no live
% visualisation, the given number of iterations). Everything the check
% changes on obj (also the fiber counts calcstats stores for the GUI) and
% the random number generator are restored afterwards.
%
% obj  - ea_disctract object.
% cvb  - from ea_permtest_cv_adapter_disctract(obj, ...).
% tol  - maximum allowed absolute difference. Default 1e-9.
%
% report - struct: .ok, .maxdiff, .NaNmatch, .Npatients, .Ncompared,
%          .explorerPred, .enginePred, .seconds (explorer, engine), .fast

if ~exist('tol', 'var') || isempty(tol)
    tol = 1e-9;
end

% the explorer's CV
props = {'kIter', 'customselection', 'useExternalModel', 'cvlivevisualize', 'nestedLOO', ...
    'adj_scaler', 'CleartuneOptim', 'stats', 'subscore'};
saved = struct;
for p = props
    saved.(p{1}) = obj.(p{1});
end
stream = rng;
cleanup = onCleanup(@() local_restore(obj, saved, stream));
obj.kIter = cvb.cv.kIter;
obj.customselection = [];
obj.useExternalModel = false;
obj.cvlivevisualize = 0;
obj.nestedLOO = false;
obj.CleartuneOptim = 0;
fprintf('Running the explorer''s own k-fold cross-validation (k = %d, %d iteration(s))...\n', cvb.cv.kfold, cvb.cv.kIter);
t0 = tic;
[~, Iexpl] = obj.kfoldcv(1);
tExpl = toc(t0);
clear cleanup

% the engine on the real scores
fprintf('Running the permutation engine''s k-fold cross-validation on the real scores...\n');
t0 = tic;
out = ea_permtest_cv_predict(cvb, ea_permtest_cv_scores(cvb));
tEng = toc(t0);

Iexpl = double(Iexpl(:));
Ieng = out.combPred(:, 1);
report = struct;
report.Npatients = numel(Ieng);
report.NaNmatch = isequal(isnan(Iexpl), isnan(Ieng));
both = ~isnan(Iexpl) & ~isnan(Ieng);
report.Ncompared = sum(both);
report.maxdiff = max([0; abs(Iexpl(both) - Ieng(both))]);
report.ok = numel(Iexpl) == numel(Ieng) && report.NaNmatch && report.maxdiff <= tol;
report.explorerPred = Iexpl;
report.enginePred = Ieng;
report.seconds = [tExpl, tEng];
report.fast = out.fast;
if report.ok
    fprintf('Check PASSED: %d patients, max |diff| %.2e (explorer %.0f s, engine %.0f s)\n', ...
        report.Ncompared, report.maxdiff, tExpl, tEng);
else
    ea_cprintf('CmdWinErrors', 'Check FAILED: %d of %d patients compared, NaN pattern equal: %d, max |diff| %.2e\n', ...
        report.Ncompared, report.Npatients, report.NaNmatch, report.maxdiff);
end


function local_restore(obj, saved, stream)
fn = fieldnames(saved);
for i = 1:numel(fn)
    obj.(fn{i}) = saved.(fn{i});
end
rng(stream);
