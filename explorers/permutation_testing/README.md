# Permutation testing (fiber filtering)

Additive add-on: nothing in the explorers' own code depends on this folder.

- `engine/`: works for any explorer tool and never touches an explorer object.
  - `ea_permtest_make_permset` / `ea_permtest_get_permset`: cohort-level shuffle set, saved to `<leadgroup dir>/permsets/`. It stores random orderings of all patients, not shuffled scores.
  - `ea_permtest_derive_shuffles`: shuffles for one patient selection. Several outcomes (subscores) are shuffled jointly, so a patient's whole row moves together.
  - `ea_permtest_bundle_shuffles`: shuffles for every analysis of an adapter bundle.
    - Subscores share one joint shuffle.
    - Each group gets the same shuffles it would get if run on its own.
  - `ea_permtest_fibers_outcome`: expands one outcome column to per-side outcome vectors (sides, mirroring), as the explorers' calcstats do.
  - `ea_permtest_fibers_run`: computes the real and the shuffled fiber maps for every analysis and side with the explorer's own test, and writes them to a new run folder in batches.
  - `ea_permtest_fastcorr_prepare` / `ea_permtest_fastcorr` / `ea_permtest_midranks`: fast path for Spearman and Pearson. The fibers are ranked once, and all shuffles of a batch are computed with one matrix product (about 180 times faster). Spearman R and p are identical to `corr` bit for bit; Pearson agrees up to rounding (~1e-15). Before use, every run checks the fast path against the test function on the real scores and the first shuffles. Other tests use the test function.
  - `ea_permtest_setup`: named analysis folders (see below).
  - `ea_permtest_cv_predict` / `ea_permtest_cv_scores` / `ea_permtest_cv_run`: prediction permutation test of the k-fold cross-validation (see below).
  - `ea_permtest_parpool`: worker count adapted to the machine (default cores − 1, at most 4; serial without the Parallel Computing Toolbox). Results do not depend on it.
  - `ea_permtest_readme_append`: the run folder's README log.
  - `ea_permtest_stats`: **shared statistics for any modality (fibers or voxels)**. It computes:
    - the Eisenstein omnibus test (Q and Q+);
    - the max-statistic FWE correction, single-step and step-down (Westfall–Young), per side, for both sides, and across all analyses;
    - the optional GPD tail fit.

    It reads the shuffled maps in chunks, so memory does not grow with map size.
  - `ea_permtest_gpd_pval`: GPD tail fit (copy of the sweetspot version).
  - `ea_permtest_fibers_stats`: one call for a fiber run folder (reader → statistics → export).
- `readers/`: describe a finished run as a "map set" for the statistics. Currently `ea_permtest_read_fiberrun`; the sweetspot reader will follow.
- `exporters/`: write the significant units. Currently `ea_permtest_export_fibers`, in the explorer's `fibcell`/`vals` format, so `ea_discfibers2trk` can read it.
- `plotting/`: null-distribution panels in the sweetspot publication style.
- `gui/ea_permtest_gui_disctract`: the **Permutation testing** menu of the fiber filtering explorer:
  - Check settings / preview
  - Run permutation test (asks for an analysis name)
  - Compute statistics
  - Prediction permutation test (k-fold CV): k, iterations and post-hoc correction come from the explorer's CV panel
  - Show / Hide significant fibers

  `ea_discfiberexplorer.mlapp` only calls `ea_permtest_gui_disctract('install', app)` from its startup function, and only if this folder is on the path.
- `adapters/`: one per explorer tool. Each extracts a neutral input bundle and checks that it reproduces the tool's own fiber values.
  - `ea_permtest_fibers_disctract`: one call for the fiber filtering explorer (adapter → check → shuffle set → run).
  - `ea_permtest_cv_disctract`: one call for the prediction permutation test (CV adapter → check against the explorer's own `kfoldcv` → shuffle set → run).
  - `ea_permtest_cv_adapter_disctract` / `ea_permtest_cv_check_disctract`: k-fold settings, fold splits and fiber values of the explorer; the check requires the engine's real predictions to match the explorer's.
  - `ea_permtest_adapter_disctract` / `ea_permtest_check_disctract`: fiber filtering explorer.
  - Supported modes: Single Tract Analysis, Split & Color By Subscore, Split & Color By Group.
  - Refused: covariates, PCA mode, hemiscores, one-sample and descriptive tests.
  - Unified explorer: to follow.

Results are saved below the lead group folder, in a folder named by the user: `<leadgroup dir>/fiberfiltering/permutation/<name>/`. One name is one analysis setup (patients, scores, fiber settings, statistic, shuffles). Its runs sit side by side:
- `README.md` and `setup.mat`: the settings of the setup. A run is only added when its settings match; otherwise a new name is needed. Existing runs are never overwritten.
- `fibers/`: fiber run (omnibus, max statistic), containing:
  - `README.md`
  - `run_info.mat`
  - `fibers.mat`
  - one `results_<subscore or group>.mat` per analysis, holding the real and shuffled R and p per side
  - `stats/` and `exports/` after "Compute statistics"
- `prediction_k<k>[_iter<n>]/`: prediction run, containing `README.md`, `run_info.mat`, `results.mat` and null-distribution figures.

Runs from before named analyses stay in `<leadgroup dir>/fiberfiltering/permutation/<ID>/` and are still listed in the GUI.

Prediction permutation test (k-fold CV): the explorer's whole k-fold pipeline is repeated for every shuffle, with the explorer's fold splits. Each fold does: coverage filter, fiber test, optional significance masking, display threshold, optional model normalization, and the prediction of the test patients. After the folds come the optional fit to scores, the subscore weights and the mean over iterations. Reported with Spearman and Pearson R (one-sided p):
- every model against its own score (uncorrected for the number of subscores);
- optionally, the combined prediction against a chosen clinical variable, which is shuffled together with the subscores.

Not supported: VTA models, 'Split & Color By Group', fiber adjacency smoothing.

```matlab
fibfiltPath = '<leadgroup dir>/fiberfiltering/<ID>.fibfilt';
S = load(fibfiltPath, '-mat'); obj = S.tractset;
run = ea_permtest_fibers_disctract(obj, struct('fibfiltPath', fibfiltPath, 'name', 'OCD_N32_pooled'));   % options: stratifyByGroup, Nperm, fastPath
cvrun = ea_permtest_cv_disctract(obj, struct('fibfiltPath', fibfiltPath, 'name', 'OCD_N32_pooled', ...
    'kIter', 1, 'targetLabel', ''));                                              % options: posthoccorrect, stratifyByGroup, Nperm
```

Statistics once the run has finished. They are written to `stats/<settings>_<time>/` and `exports/<same>/`, and summarised in the run's README:

```matlab
[stats, files] = ea_permtest_fibers_stats(run.rundir);   % defaults: positive direction, alpha 0.05, all families, Q and Q+
```

The same steps one by one:

```matlab
bundle = ea_permtest_adapter_disctract(obj, struct('fibfiltPath', fibfiltPath));
report = ea_permtest_check_disctract(obj, bundle);
permset = ea_permtest_get_permset(bundle.groupdir, bundle.patientIDs, 1000, bundle.rngseed);
shuffles = ea_permtest_bundle_shuffles(permset, bundle);   % one entry per subscore / group / single analysis
run = ea_permtest_fibers_run(bundle, permset, struct('name', 'OCD_N32_pooled', 'checkReport', report));
```
