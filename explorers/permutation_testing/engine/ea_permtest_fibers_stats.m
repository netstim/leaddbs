function [stats, exportFiles] = ea_permtest_fibers_stats(rundir, opts)
% Statistics for a finished fiber permutation run (ea_permtest_fibers_run):
% Eisenstein omnibus tests and FWE-corrected max statistic (single-step
% and step-down), computed by the shared engine ea_permtest_stats, plus
% export of the significant fibers for 3D visualisation.
%
% Fiber defaults (all can be changed through opts, see ea_permtest_stats):
%   direction 'positive' (improvement-related fibers), alpha 0.05,
%   families perside + bothsides + allanalyses, omnibus Q and Qpos per
%   side and for both sides, step-down on, GPD tail fit off.
%
% Outputs (never overwriting earlier statistics):
%   <rundir>/stats/<direction>_alpha<alpha>_<time>/   stats.mat + figures
%   <rundir>/exports/<same name>/                     significant fibers
%   <rundir>/README.md                                summary appended
%
% rundir - fiber permutation run folder (run.rundir).
% opts   - optional struct, passed to ea_permtest_stats.

if ~exist('opts', 'var') || isempty(opts)
    opts = struct;
end
mapset = ea_permtest_read_fiberrun(rundir);
stats = ea_permtest_stats(mapset, opts);
[~, name, ext] = fileparts(stats.settings.outdir);
name = [name, ext]; % the folder name contains the alpha level, e.g. ..._alpha0.05_<time>
exportFiles = ea_permtest_export_fibers(stats, mapset, fullfile(mapset.exportsRoot, name));
