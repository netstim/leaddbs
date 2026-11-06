function mapset = ea_permtest_read_fiberrun(rundir)
% Describe a finished fiber permutation run folder (ea_permtest_fibers_run)
% as a "map set" for the shared statistics engine (ea_permtest_stats).
% Only the small per-unit vectors are loaded; the large shuffled matrices
% stay on disk and are read in chunks by the engine.
%
% A map is one analysis (single tract / subscore / group) on one side; its
% units are that side's analysed fibers.
%
% mapset fields:
%   modality    'fibers'
%   source      rundir
%   analyses(a) .label, .fiberset
%   maps(i)     .analysis, .side, .sideLabel, .label, .Nunits, .Nperm,
%               .Remp, .pemp (Nunits x 1), .unitIdx (explorer fiber rows),
%               .file, .Rvar, .pvar (variables holding the Nperm x Nunits
%               shuffled statistics / p-values in .file)
%   readmeDir   folder whose README.md documents the run
%   statsRoot, exportsRoot   where statistics / exports are written
%   meta        .fibersFile (fiber coordinates), .settings (run settings)

if ~isfolder(rundir)
    ea_error(['Run folder not found:\n', strrep(strrep(rundir, '\', '\\'), '%', '%%')]);
end
infoFile = fullfile(rundir, 'run_info.mat');
if ~isfile(infoFile)
    ea_error('No run_info.mat in this folder -- not a fiber permutation run folder.');
end
ri = load(infoFile);
if ~isfield(ri, 'seconds')
    ea_error(['This permutation run has not finished (or was interrupted): ', ...
        'run_info.mat has no timing entry yet. Statistics need a complete run.']);
end
settings = ri.settings;

% Result files, matched to the analyses by their stored label
files = dir(fullfile(rundir, 'results_*.mat'));
labels = cell(1, numel(files));
for f = 1:numel(files)
    labels{f} = matfile(fullfile(rundir, files(f).name)).label;
end

sideLabels = {'right', 'left'};
mapset = struct;
mapset.modality = 'fibers';
mapset.source = rundir;
mapset.analyses = struct('label', {}, 'fiberset', {});
mapset.maps = struct('analysis', {}, 'side', {}, 'sideLabel', {}, 'label', {}, 'Nunits', {}, ...
    'Nperm', {}, 'Remp', {}, 'pemp', {}, 'unitIdx', {}, 'file', {}, 'Rvar', {}, 'pvar', {});
for a = 1:numel(settings.analyses)
    label = settings.analyses(a).label;
    f = find(strcmp(labels, label), 1);
    if isempty(f)
        ea_error(sprintf('No results file for analysis ''%s'' in this run folder.', label));
    end
    file = fullfile(rundir, files(f).name);
    m = matfile(file);
    mapset.analyses(a).label = label;
    mapset.analyses(a).fiberset = m.fiberset;
    for side = 1:settings.Nsides
        Rvar = sprintf('Rperm_side%d', side);
        sz = size(m, Rvar);
        if sz(2) == 0
            continue % no analysed fibers on this side
        end
        lastRow = m.(Rvar)(sz(1), :);
        if all(isnan(lastRow))
            ea_error(sprintf('Results for ''%s'' (side %d) are incomplete -- the run did not finish.', label, side));
        end
        i = numel(mapset.maps) + 1;
        mapset.maps(i).analysis = a;
        mapset.maps(i).side = side;
        mapset.maps(i).sideLabel = sideLabels{min(side, 2)};
        mapset.maps(i).label = sprintf('%s, %s', label, mapset.maps(i).sideLabel);
        mapset.maps(i).Nunits = sz(2);
        mapset.maps(i).Nperm = sz(1);
        mapset.maps(i).Remp = m.(sprintf('Remp_side%d', side));
        mapset.maps(i).pemp = m.(sprintf('pemp_side%d', side));
        mapset.maps(i).unitIdx = m.(sprintf('fibIdx_side%d', side));
        mapset.maps(i).file = file;
        mapset.maps(i).Rvar = Rvar;
        mapset.maps(i).pvar = sprintf('pperm_side%d', side);
    end
end

mapset.readmeDir = rundir;
mapset.statsRoot = fullfile(rundir, 'stats');
mapset.exportsRoot = fullfile(rundir, 'exports');
mapset.meta = struct('fibersFile', fullfile(rundir, 'fibers.mat'), 'settings', settings);
