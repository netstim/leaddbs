function permresults = ea_unified_perm_sweetspot(obj, Nperm, corrType, stratifyByGroup, maxWorkers)
% Voxelwise permutation test of the sweetspotmapping correlation map.
% Permutes obj.responsevar (optionally group-restricted, NaN-excluded,
% tied across sides/mirrors since permutation happens before the L/R &
% mirror expansion in ea_unified_perm_sweetspot_calcstats), rebuilds the
% map for each permutation, and saves the resulting empirical + permuted
% voxelwise R/p maps together with the settings used, in their own
% dedicated output folder.
%
% obj             - an ea_unifiedmapping object with sweetspotmapping
%                    results already computed (obj.results.sweetspotmapping).
% Nperm           - number of permutations. Default 1000.
% corrType        - correlation type passed to ea_corr. Default
%                    obj.statsettings.stattest (the same test the live map
%                    already uses), so permtest matches obj.draw() by
%                    default unless you deliberately override it.
% stratifyByGroup - if true, shuffles obj.responsevar only within
%                    obj.M.patient.group buckets. Default false.
% maxWorkers      - cap on parpool workers for the parfor loop below (each
%                    worker holds a full copy of obj, including
%                    obj.results.sweetspotmapping.efield). Default 3.
%
% None of these settings are stored on obj -- they're plain arguments here
% (and app-level GUI properties, not model properties, when driven from the
% .mlapp), so this feature adds nothing to the size of a saved .explorer
% file. All settings actually used for this run, plus the results
% themselves, are saved together in one self-contained output folder (see
% outdir below) -- nothing else on obj is touched or grown.

if ~exist('Nperm','var') || isempty(Nperm)
    Nperm = 1000;
end
if ~exist('corrType','var') || isempty(corrType)
    corrType = obj.statsettings.stattest;
end
if ~exist('stratifyByGroup','var') || isempty(stratifyByGroup)
    stratifyByGroup = false;
end
if ~exist('maxWorkers','var') || isempty(maxWorkers)
    maxWorkers = 3;
end

% This feature only supports sweetspotmapping today. Checked explicitly
% (not inferred from obj.drawTool) because obj.drawTool reflects whatever
% is currently displayed/active in the GUI, not necessarily what the
% caller means to test -- an object can have sweetspot AND fiberfiltering/
% networkmapping results computed at once, with drawTool currently pointed
% at one of the others. ea_unified_perm_sweetspot_calcstats.m unconditionally
% reads obj.results.sweetspotmapping.efield regardless of drawTool, so
% without this check a caller could silently get a valid-looking sweetspot
% permutation test while believing they tested a different tool entirely.
if ~isfield(obj.results, 'sweetspotmapping')
    ea_error(['permtest currently only supports sweetspotmapping (fiberfiltering/networkmapping ', ...
        'are not yet implemented). obj.results.sweetspotmapping was not found -- compute ', ...
        'sweetspot results for this object first.']);
end

if ~strcmp(obj.statsettings.statfamily, 'Correlations') || ...
        ~ismember(obj.statsettings.stimulationmodel, {'Electric Field','Sigmoid Field'})
    ea_error(['permtest is currently only implemented for statsettings.statfamily = ''Correlations'' ', ...
        'with statsettings.stimulationmodel in {''Electric Field'',''Sigmoid Field''}.']);
end

if size(obj.responsevar, 2) > 1
    ea_error('Hemiscore responsevar (2 columns) is not yet supported for permutation testing.');
end

patsel = obj.patientselection;

if stratifyByGroup
    groupvec = obj.M.patient.group;
else
    groupvec = [];
end

[Iperm, PermIdx] = ea_shuffle_grouped(obj.responsevar, Nperm, patsel, groupvec, obj.rngseed);

% Descriptive folder name -- encodes the options that most affect the
% result, plus a timestamp so repeated runs never silently collide. Every
% file for this run (results, its slim companion once created, exported
% niftis, plots) lives inside this one folder -- nothing shared/flat
% across runs, nothing written outside it.
if stratifyByGroup
    stratTag = 'stratified';
else
    stratTag = 'pooled';
end
if obj.mirrorsides
    mirrorTag = 'mirrored';
else
    mirrorTag = 'nonmirrored';
end
runName = sprintf('%s_permtest_N%d_%s_%s_%s_%s', obj.ID, Nperm, stratTag, mirrorTag, corrType, datestr(now, 'yyyymmdd_HHMMSS'));
outdir = fullfile(fileparts(obj.leadgroup), 'UnifiedMappingExplorer', 'permutation', runName, filesep);
suffix = 1;
while exist(outdir, 'dir') % timestamp makes this vanishingly rare, but never silently overwrite/merge into an existing run's folder
    suffix = suffix + 1;
    outdir = fullfile(fileparts(obj.leadgroup), 'UnifiedMappingExplorer', 'permutation', sprintf('%s_%d', runName, suffix), filesep);
end
ea_mkdir(outdir);

fprintf('Calculating empirical (unpermuted) sweetspot map...\n');
[Remp, pemp, gvalFixed, gpatselFixed, thisvalsFixed, nanidxFixed] = ea_unified_perm_sweetspot_calcstats(obj, patsel, obj.responsevar, true); % skipsigthresh=true: always store raw, unmasked values
% gvalFixed/gpatselFixed (coverage-masked efield data & patient selection)
% and thisvalsFixed/nanidxFixed (the patient-sliced, NaN-filtered
% correlation input derived from them) don't depend on I/Iperm -- reused
% for every permutation below instead of being recomputed (which forces
% full-matrix copies) on every one of the Nperm calls.

fprintf('Checking empirical map against the live explorer calcstats path...\n');
ea_unified_perm_sweetspot_consistency_check(obj, patsel, Remp);

% Export the empirical (raw, unmasked) R-map to nifti by default, so it can
% be sanity-checked against any map generated the usual way (e.g. obj.draw()).
ea_unified_perm_vals2nii(obj.results.sweetspotmapping.space, Remp, outdir, 'Remp');

ea_unified_perm_cap_parpool(maxWorkers);

Rrow = cell(1, Nperm);
prow = cell(1, Nperm);

progress = 0;
dq = parallel.pool.DataQueue;
afterEach(dq, @(~) reportProgress());

fprintf('Running %d permutations...\n', Nperm);
parfor p = 1:Nperm
    [v, pv] = ea_unified_perm_sweetspot_calcstats(obj, patsel, Iperm(:,p), true, gvalFixed, gpatselFixed, thisvalsFixed, nanidxFixed);
    Rrow{p} = v;
    prow{p} = pv;
    send(dq, 1);
end

% Rperm/pperm are assembled in memory here (double precision, one side at
% a time is not needed for this step -- this reassembly was never the
% documented memory-crash point; that was always downstream code
% blanket-loading the SAVED file back in). What actually goes to disk
% below is single-precision and flat (not nested in a cell), directly in
% the same shape downstream readers (ea_unified_perm_sweetspot_predict.m,
% ea_unified_perm_sweetspot_stats.m) need for chunked matfile reads -- so
% there is exactly one copy of this data on disk, at half the size double
% precision would take, with no separate conversion step required before
% first use. Remp/pemp (small, Nvoxels x 1 per side) stay double precision
% and nested in a friendly {side} cell -- no memory pressure motivates
% downcasting them, and ea_unified_perm_sweetspot_consistency_check.m's
% tight tolerance against the live explorer calcstats path depends on
% Remp keeping full double precision.
Rperm = cell(size(Remp));
pperm = cell(size(Remp));
for side = 1:numel(Remp)
    Rperm{side} = nan(Nperm, numel(Remp{side}));
    pperm{side} = nan(Nperm, numel(Remp{side}));
    for p = 1:Nperm
        Rperm{side}(p,:) = Rrow{p}{side}';
        pperm{side}(p,:) = prow{p}{side}';
    end
end

permresults.Remp = Remp;
permresults.pemp = pemp;
permresults.nsides = numel(Remp);
for side = 1:numel(Remp)
    permresults.(sprintf('Rperm_side%d', side)) = single(Rperm{side});
    permresults.(sprintf('pperm_side%d', side)) = single(pperm{side});
end
permresults.PermIdx = PermIdx;
permresults.Nperm = Nperm;
permresults.corrType = corrType;
permresults.space = obj.results.sweetspotmapping.space; % training grid geometry, needed for out-of-sample reslicing later
permresults.ID = obj.ID;
permresults.leadgroup = obj.leadgroup;
permresults.patientselection = patsel;
permresults.mirrorsides = obj.mirrorsides;
permresults.stratifyPermutationsByGroup = stratifyByGroup;

% Full record of every setting this run actually used -- settings and
% results kept together in one self-contained file/folder, not split
% across the model object and disk.
permresults.settings.Nperm = Nperm;
permresults.settings.corrType = corrType;
permresults.settings.stratifyByGroup = stratifyByGroup;
permresults.settings.maxWorkers = maxWorkers;
permresults.settings.rngseed = obj.rngseed;
permresults.settings.mirrorsides = obj.mirrorsides;
permresults.settings.statfamily = obj.statsettings.statfamily;
permresults.settings.stimulationmodel = obj.statsettings.stimulationmodel;
permresults.settings.efieldthreshold_spot = obj.statsettings.efieldthreshold_spot;
permresults.settings.connthreshold = obj.statsettings.connthreshold;
permresults.settings.nanthreshold = obj.statsettings.nanthreshold;
permresults.settings.showsignificantonly = obj.showsignificantonly;
permresults.settings.multcompstrategy = obj.multcompstrategy;
permresults.settings.alphalevel = obj.alphalevel;
permresults.settings.hasCovars = ~isempty(obj.covars);

outfile = fullfile(outdir, 'results.mat');
permresults.savedfile = outfile;
save(outfile, '-struct', 'permresults', '-v7.3');
fprintf('Saved permutation results to %s\n', outfile);

readmeBody = sprintf([ ...
    'Voxelwise permutation test of the sweetspotmapping correlation map for %s.\n\n', ...
    'Settings used for this run:\n', ...
    '  Nperm                = %d\n', ...
    '  corrType              = %s\n', ...
    '  stratifyByGroup       = %d\n', ...
    '  mirrorsides           = %d\n', ...
    '  statfamily            = %s\n', ...
    '  stimulationmodel      = %s\n', ...
    '  efieldthreshold_spot  = %g\n', ...
    '  connthreshold         = %g%%\n', ...
    '  rngseed               = %s\n\n', ...
    'Files in this folder:\n', ...
    '  results.mat  - Remp/pemp (empirical R/p per voxel, double precision) and\n', ...
    '                 Rperm_side*/pperm_side* (Nperm x Nvoxels null distributions,\n', ...
    '                 single precision, one flat variable per side -- read these via\n', ...
    '                 matfile() row/column-chunked access, do not load() them whole),\n', ...
    '                 plus PermIdx (the exact patient-shuffle used per permutation)\n', ...
    '                 and the settings block above.\n', ...
    '  Remp_r.nii / Remp_l.nii - the empirical (real, unpermuted) correlation map,\n', ...
    '                 UNMASKED (raw R everywhere, no significance thresholding applied\n', ...
    '                 yet -- that happens in a separate post-hoc step, see\n', ...
    '                 ea_unified_perm_sweetspot_stats, which adds its own files/section\n', ...
    '                 to this same folder/README once run).\n'], ...
    obj.ID, Nperm, corrType, stratifyByGroup, obj.mirrorsides, obj.statsettings.statfamily, ...
    obj.statsettings.stimulationmodel, obj.statsettings.efieldthreshold_spot, obj.statsettings.connthreshold, obj.rngseed);
ea_unified_perm_readme_append(outdir, 'Voxelwise permutation test (ea_unified_perm_sweetspot)', readmeBody);

    function reportProgress()
        % Called on the client (not inside the workers) each time a
        % worker finishes one permutation, via the DataQueue above.
        progress = progress + 1;
        step = max(1, round(Nperm/20)); % ~5% increments
        if mod(progress, step) == 0 || progress == Nperm
            fprintf('Permutation progress: %d/%d (%.0f%%)\n', progress, Nperm, 100*progress/Nperm);
        end
    end
end
